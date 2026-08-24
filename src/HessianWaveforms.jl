using LinearAlgebra

#
# Low-rank Hessian analysis, part 2: eliminating nuisance parameters, the Hessian of the
# ROBUSTNESS cost, and turning eigenvectors into per-quadrature waveforms.
#
# `HessianAnalysis.jl` answers "how does the fidelity degrade under a perturbation of the
# controls". That is the object a closed-loop calibration searches in (arXiv:2606.05060).
# Two things are missing before it can be used that way:
#
#   1. A calibration never holds the free single-qubit phase fixed -- it recalibrates it.
#      The relevant curvature is therefore that of `max_{x_add} F`, which is the SCHUR
#      COMPLEMENT of the full Hessian on the additional parameters, not its main block.
#
#   2. Cancelling a distortion restores the fidelity at the nominal operating point but
#      says nothing about robustness. Grading the double Taylor series of F in (error eps,
#      control perturbation du) by (order in eps, order in du):
#
#          (0,2)  fidelity Hessian H                                    <- HessianAnalysis
#          (1,1)  ZERO: F <= 1 is attained, so the joint Hessian is negative semidefinite,
#                 and a PSD matrix with a zero diagonal entry (robustness => d2F/deps2 = 0)
#                 has that whole row and column zero.
#          (2,1)  ZERO: writing the robustness cost as A(u) = 2||E(u)||^2 with E the
#                 first-order error operator, grad A = 4 Re<E, dE/du> vanishes at E = 0.
#          (2,2)  K = d^2A/du^2 = 4 G'G,  G = dE/du.                    <- THIS FILE
#
#      So K -- the mixed fourth derivative d^4F/(deps^2 du_i du_j) -- is the LEADING term
#      governing how much robustness a correction costs, not a subleading one.
#
#      K IS A GRAM MATRIX ONLY AT A ROBUST POINT. The full second derivative is
#      d^2A = 4[Re<E, d^2E/du^2> + G'G]; the first term drops out precisely when E = 0, and
#      only then is K PSD and of the same low rank as the error-operator space. Evaluated at
#      a NON-robust pulse -- the time-optimal gate, say -- K is an ordinary symmetric matrix
#      with negative eigenvalues, and that is correct rather than a numerical failure.
#
# K IS NOT COMPUTED AS A FOURTH DERIVATIVE. Averaging the infidelity over +-sigma of the
# error kills the odd terms, so
#
#     R_sigma(u) = (1/2)[(1 - F(u,+sigma)) + (1 - F(u,-sigma))]
#                = (1 - F(u,0)) + (1/2) sigma^2 A(u) + O(sigma^4)
#
# and therefore  d^2R_sigma/du^2 = H + (1/2) sigma^2 K + O(sigma^4). Both Hessians on the
# right come from `calculate_fidelity_hessian` evaluated at a STATIC error offset, so no
# new derivative code is involved and the existing validation carries over. It is also
# what an experiment actually measures: fidelity at a deliberately detuned operating
# point. `calculate_robustness_hessian` returns K = 2(R'' - H)/sigma^2.
#
# ORDER OF OPERATIONS MATTERS. The Schur complement is the Hessian of the partially
# maximized function only where the gradient w.r.t. the eliminated parameters vanishes,
# and the optimal x_add MOVES with sigma. So x_add is re-optimized at each sigma and the
# Schur complement is taken BEFORE the +-sigma average, never after.
#

export schur_eliminate, optimize_additional_parameters, shift_error_source
export calculate_robustness_hessian, principal_waveforms, hessian_health
export apply_control_constraint

"""
    schur_eliminate(H::AbstractMatrix, idx) -> Matrix

Curvature that survives after the parameters in `idx` are re-optimized.

For `H = -∂²F/∂x∂x` partitioned into kept (`k`) and eliminated (`e`) blocks, this returns
``H_{kk} - H_{ke} H_{ee}^{-1} H_{ek}``, which is ``-∂²/∂x_k²`` of ``\\max_{x_e} F``.

This is only the Hessian of the maximized function if ``∂F/∂x_e = 0`` at the expansion
point; the caller is responsible for that (see [`optimize_additional_parameters`]).

`idx` may be empty, in which case `H` is returned unchanged.

# Arguments
- `H`: symmetric matrix
- `idx`: indices to eliminate

# Errors
Throws if the eliminated block is numerically singular, which means one of those
parameters is a flat direction and cannot be "optimized away" -- silently pseudo-inverting
there would return a confident wrong curvature.
"""
function schur_eliminate(H::AbstractMatrix, idx)
    n = size(H, 1)
    size(H, 2) == n || throw(ArgumentError("schur_eliminate: H must be square"))
    idx = collect(idx)
    isempty(idx) && return Matrix(H)
    all(1 .<= idx .<= n) || throw(ArgumentError("schur_eliminate: idx out of range"))
    keep = setdiff(1:n, idx)
    A = H[keep, keep]
    B = H[keep, idx]
    D = Symmetric(Matrix(H[idx, idx]))
    rc = rcond_sym(D)
    rc > 1e-10 || error("schur_eliminate: eliminated block is singular (rcond = $rc). " *
                        "One of the parameters $(idx) is a flat direction; it cannot be " *
                        "re-optimized away.")
    return Matrix(Symmetric(A - B * (D \ transpose(B))))
end

"""Reciprocal condition number of a symmetric matrix, via its eigenvalues."""
function rcond_sym(D::Symmetric)
    λ = abs.(eigvals(D))
    m = maximum(λ)
    m == 0 && return 0.0
    return minimum(λ) / m
end

"""
    optimize_additional_parameters(fp, x; kwargs...) -> Vector

Newton refinement of the additional parameters (the trailing `nb_additional_param` entries
of `x`) at fixed control waveform, using the analytic fidelity gradient.

This is a LOCAL refiner: the fidelity of a CZ target is a trigonometric polynomial in the
free single-qubit phase and generally has several stationary points, so `x` must already
be near the right one. It is meant for tracking how the optimum MOVES when the operating
point is shifted (see [`calculate_robustness_hessian`]), not for finding it cold.

# Arguments
- `fp::FidelityRobustGRAPEProblem`
- `x`: full parameter vector, whose additional entries are used as the seed

# Keywords
- `iterations`: maximum Newton steps (default 50)
- `gtol`: stop once the additional-parameter gradient is below this (default 1e-9). The
  floor is not adjustable downward in practice: `calculate_fidelity_and_derivatives`
  forward-differences the TARGET unitary w.r.t. `x_add` with step `UnitaryRobustGRAPEProblem.ϵ`
  (1e-8 by default), so `dF/dx_add` is itself only accurate to ~1e-9.
- `steptol`: also stop once the Newton step is below this, which is what actually
  terminates the iteration once the gradient has hit its own noise floor
- `h`: finite-difference step for the Newton matrix (default 1e-6)
- `verbose`: print each step

# Returns
A copy of `x` with the additional parameters refined.
"""
function optimize_additional_parameters(fp::FidelityRobustGRAPEProblem, x::Vector{<:Real};
                                        iterations::Int=50, gtol::Real=1e-9,
                                        steptol::Real=1e-12, h::Real=1e-6,
                                        verbose::Bool=false)
    nb = fp.unitary_problem.nb_additional_param
    nb == 0 && return copy(x)
    N = length(x)
    idx = (N - nb + 1):N
    xc = collect(float.(x))
    gadd(v) = calculate_fidelity_and_derivatives(fp, v)[2][idx]

    for it in 1:iterations
        g = gadd(xc)
        if norm(g) < gtol
            verbose && println("  x_add converged in $(it - 1) steps, |g| = $(norm(g))")
            return xc
        end
        M = zeros(nb, nb)
        for j in 1:nb
            xp = copy(xc); xp[idx[j]] += h
            xm = copy(xc); xm[idx[j]] -= h
            M[:, j] = (gadd(xp) - gadd(xm)) ./ (2h)
        end
        M = (M + transpose(M)) ./ 2
        step = -(M \ g)
        # A Newton step toward a SADDLE is possible if the curvature block is indefinite;
        # fall back to gradient ascent there rather than walking off the maximum.
        if !all(eigvals(Symmetric(M)) .< 0)
            step = g ./ max(norm(g), eps())
        end
        xc[idx] .+= step
        verbose && @info "x_add Newton" it norm(g) step
        norm(step) < steptol && return xc
    end
    @warn "optimize_additional_parameters: did not converge in $iterations steps" resid=norm(gadd(xc))
    return xc
end

"""
    shift_error_source(fp, k, σ) -> FidelityRobustGRAPEProblem

The same problem evaluated at a STATIC offset `σ` of error source `k`: the offset is
folded into `H0` and the error sources are dropped.

Exact, not perturbative, PROVIDED `Herror(nt, x, x_add, ε)` means "the change in the
Hamiltonian at error amplitude ε" -- the documented convention for [`ErrorSource`]. That
is checked at construction by evaluating `Herror(..., 0)`, which must vanish.
"""
function shift_error_source(fp::FidelityRobustGRAPEProblem, k::Int, σ::Real;
                            keep_sources::Bool=false)
    up = fp.unitary_problem
    1 <= k <= length(up.error_sources) ||
        throw(ArgumentError("shift_error_source: no error source $k"))
    Hk = up.error_sources[k].Herror
    H0old = up.H0
    H0new = HamiltonianFunctionWrapper(
        (nt, x, xa) -> H0old(nt, x, xa) .+ Hk(nt, x, xa, σ))
    up2 = UnitaryRobustGRAPEProblem(
        t0=up.t0, ntimes=up.ntimes, ndim=up.ndim, H0=H0new,
        nb_additional_param=up.nb_additional_param,
        error_sources=(keep_sources ? up.error_sources : ErrorSource[]),
        ϵ=up.ϵ, ϵ2=up.ϵ2)
    return FidelityRobustGRAPEProblem(unitary_problem=up2, projector=fp.projector,
                                      target_unitary=fp.target_unitary)
end

"""
    calculate_robustness_hessian(fp, x; kwargs...) -> NamedTuple

Curvature of the ROBUSTNESS cost with respect to the control waveform, alongside the
ordinary fidelity Hessian, both with the additional parameters eliminated.

Let ``A(u) = -∂²F/∂ε²`` be the first-order sensitivity to error source `source` and
``\\tilde F = \\max_{x_{add}} F``. Averaging the infidelity over ``\\pm σ`` kills the odd
terms of the expansion in ε, so

```math
R_σ(u) = \\tfrac12\\left[(1-\\tilde F(u,σ)) + (1-\\tilde F(u,-σ))\\right]
       = (1-\\tilde F(u,0)) + \\tfrac12 σ^2 A(u) + O(σ^4)
```

and this function returns ``K = 2(R_σ'' - H)/σ^2 = ∂^2A/∂u^2 + O(σ^2)``, the mixed fourth
derivative ``∂^4F/(∂ε^2 ∂u_i ∂u_j)``.

`K` is positive semidefinite and low rank ONLY where the gate is actually robust
(``E = 0``), which is where it is meant to be used; at a non-robust pulse it is an ordinary
symmetric matrix and negative eigenvalues are physical, not numerical.

`σ` is a numerical parameter, NOT a physical choice: it cancels out of `K` to `O(σ²)`, and
in particular the eigenVECTORS of `K` do not depend on it. Verify by halving it -- `K`
must be unchanged. Too small and the difference of two nearly equal Hessians is dominated
by the finite-difference floor of each; too large and the `O(σ²)` truncation shows. The
default sits near the optimum of that trade-off for well-conditioned problems.

# Arguments
- `fp::FidelityRobustGRAPEProblem`: must carry at least one error source
- `x`: parameter vector at the design point (`F ≈ 1`)

# Keywords
- `source`: which error source to expand in (default 1)
- `σ`: static offset used for the difference (default 0.05)
- `ε_hess`: passed through to [`calculate_fidelity_hessian`]
- `reoptimize_add`: re-optimize the additional parameters at each ±σ (default `true`).
  Required for the Schur complement to mean anything -- see the note on ordering below.
- `constraint`: optional `N × m` matrix tying control parameters together before the Schur
  step, for problems whose physical knobs are fewer than the propagation's parameters --
  a composite gate that plays ONE waveform in several segments, say. See
  [`apply_control_constraint`] for the required shape.
- `verbose`

# Returns
NamedTuple with
- `K`: the robustness-cost Hessian, additional parameters eliminated
- `fidelity_hessian`: `H` at σ = 0, additional parameters eliminated
- `robustness_cost_hessian`: ``R_σ''``, i.e. `H + σ²K/2`
- `sigma`, `x_add`: the σ used and the re-optimized additional parameters at (0, +σ, −σ)
- `A0`: the sensitivity ``-∂²\tilde F/∂ε²`` at the design point, by the same ±σ average.
  This is the CORRECTED sensitivity (x_add re-optimized), so for a pseudo-robust gate it is
  the near-zero number, not the raw one. `K` is a Gram matrix only where `A0 ≈ 0` -- check
  it before reading any structure into `K`.
- `infidelity`: `1 - F̃` at (0, +σ, −σ), for sanity

# Note on ordering
The Schur complement is the curvature of the maximized function only where the gradient
w.r.t. the eliminated parameters vanishes, and the optimum moves with σ. The additional
parameters are therefore re-optimized at each σ and eliminated BEFORE the ±σ average is
taken. Averaging first and eliminating afterwards is wrong and does not announce itself.

A `constraint` is applied before the Schur complement for the same reason: the physical
function is `F(w, x_add) = F(C·(w, x_add))`, whose Hessian is `Cᵀ H C`; only then does
eliminating `x_add` mean "re-optimize the single-qubit phase at fixed physical waveform".
"""
function calculate_robustness_hessian(fp::FidelityRobustGRAPEProblem, x::Vector{<:Real};
                                      source::Int=1, σ::Real=0.05, ε_hess::Real=1e-5,
                                      reoptimize_add::Bool=true,
                                      constraint::Union{Nothing,AbstractMatrix}=nothing,
                                      verbose::Bool=false)
    up = fp.unitary_problem
    nb = up.nb_additional_param
    N = length(x)
    idx_add = (N - nb + 1):N
    ntimes = up.ntimes
    nparam = div(N - nb, ntimes)

    # The ±σ shift is only exact if Herror is the FULL change at amplitude ε, i.e. it
    # vanishes at ε = 0.  Checked, not assumed.
    let xm2 = reshape(collect(float.(x[1:end-nb])), (nparam, ntimes)),
        xa = collect(float.(x[idx_add]))
        z = up.error_sources[source].Herror(1, xm2[:, 1], xa, 0.0)
        norm(z) < 1e-12 || error("calculate_robustness_hessian: error source $source does " *
            "not vanish at ε = 0 (‖Herror(…,0)‖ = $(norm(z))); the static-offset " *
            "construction assumes Herror(ε) is the full change at amplitude ε.")
    end

    function hess_at(s)
        fps = s == 0 ? fp : shift_error_source(fp, source, s)
        xs = (reoptimize_add && nb > 0) ? optimize_additional_parameters(fps, x) : collect(float.(x))
        Fs = calculate_fidelity_and_derivatives(fps, xs)[1]
        H = calculate_fidelity_hessian(fps, xs; ε_hess=ε_hess, verbose=verbose)
        Hc, idx_red = apply_control_constraint(H, constraint, nb)
        return (schur_eliminate(Hc, idx_red), xs[idx_add], Fs)
    end

    verbose && println("  robustness Hessian: σ = 0")
    H0, a0, F0 = hess_at(0)
    verbose && println("  robustness Hessian: σ = +$σ")
    Hp, ap, Fp = hess_at(+σ)
    verbose && println("  robustness Hessian: σ = -$σ")
    Hm, am, Fm = hess_at(-σ)

    Hrob = (Hp .+ Hm) ./ 2
    K = 2 .* (Hrob .- H0) ./ σ^2
    # A0 is the CORRECTED sensitivity at the design point -- the same ±σ average, one
    # order down. It decides whether K means anything: K = 4G'G, positive semidefinite and
    # low rank, only where A0 = 0. A gate whose robustness is claimed but whose A0 is not
    # small will produce an indefinite, full-rank K, and that is a statement about the
    # gate, not about the numerics.
    A0 = -(Fp + Fm - 2F0) / σ^2
    return (K=Matrix(Symmetric(K)), fidelity_hessian=Matrix(Symmetric(H0)),
            robustness_cost_hessian=Matrix(Symmetric(Hrob)), sigma=σ,
            A0=A0, infidelity=(1 - F0, 1 - Fp, 1 - Fm), x_add=(a0, ap, am))
end

"""
    principal_waveforms(H, ntimes, nparam; kwargs...) -> NamedTuple

Eigen-decompose a control-space Hessian and reshape each eigenvector into one waveform per
control quadrature.

The parameter layout matches the rest of the package: `x` is `reshape`d to
`(nparam, ntimes)`, so entry `(np, nt)` of a returned waveform is quadrature `np` at time
step `nt`. For a two-control Rydberg problem `p = [phase, relative modulus]`, `waveforms[k]`
is a `2 × ntimes` matrix whose rows are ``δφ_k(t)`` and ``δ|Ω|_k(t)/|Ω|``.

Polar versus Cartesian is a per-timestep ROTATION of the control basis
(``δ\\tildeΩ = e^{-iφ}(δA - iΩ δφ)``), so it leaves the eigenvalues untouched and only
re-expresses the eigenvectors. Polar is preferred because both components are then
dimensionless -- radians and relative amplitude -- and directly comparable without
choosing a scale.

# Keywords
- `n_max`: how many directions to return (default 10)
- `rtol`: report the rank as the number of eigenvalues above `rtol * λ₁` (default 1e-4)
- `degeneracy_rtol`: warn when two returned eigenvalues are closer than this in relative
  terms (default 1e-3). Within a degenerate block the individual eigenvectors are only
  defined up to a rotation, so plotting them as if they were canonical is meaningless.

# Returns
NamedTuple with `eigenvalues`, `eigenvectors` (columns), `waveforms` (vector of
`nparam × ntimes` matrices), `rank`, and `spectrum` (all eigenvalues, descending).

# Sign convention
Eigenvectors are fixed so that the component of largest magnitude is positive. Without
this the sign is whatever LAPACK returned and the figure will not reproduce.
"""
function principal_waveforms(H::AbstractMatrix, ntimes::Int, nparam::Int;
                             n_max::Int=10, rtol::Real=1e-4, degeneracy_rtol::Real=1e-3)
    N = size(H, 1)
    N == nparam * ntimes ||
        throw(ArgumentError("principal_waveforms: H is $(N)×$(N) but nparam*ntimes = " *
                            "$(nparam * ntimes). Eliminate the additional parameters " *
                            "first (see schur_eliminate)."))
    eig = eigen(Symmetric(Matrix(H)))
    perm = sortperm(eig.values; rev=true)
    λ = eig.values[perm]
    V = eig.vectors[:, perm]

    for k in 1:size(V, 2)                   # sign convention
        V[:, k] .*= sign(V[argmax(abs.(@view V[:, k])), k])
    end

    r = count(λ .> rtol * λ[1])
    n = min(n_max, length(λ))
    for k in 2:n
        if abs(λ[k] - λ[k-1]) < degeneracy_rtol * abs(λ[1])
            @warn "principal_waveforms: eigenvalues $(k-1) and $k are nearly degenerate " *
                  "($(λ[k-1]) vs $(λ[k])); their eigenvectors are only defined up to a " *
                  "rotation within that block."
        end
    end

    waveforms = [reshape(V[:, k], (nparam, ntimes)) for k in 1:n]
    return (eigenvalues=λ[1:n], eigenvectors=V[:, 1:n], waveforms=waveforms,
            rank=r, spectrum=λ)
end

"""
    hessian_health(fp, x, H; rtol=1e-6) -> NamedTuple

Cheap diagnostics for a Hessian that is supposed to sit at a fidelity maximum.

`calculate_fidelity_hessian` finite-differences `exp(-i dt H0)` and divides by `ε_hess²`.
When `‖dt·H0‖` is large -- a Rydberg problem with an explicit large blockade reaches
`B·dt ~ 5e4` -- roundoff is amplified by `1/ε_hess²` and the routine returns a matrix that
is not even positive semidefinite, without raising anything. This checks the two things
that must hold at a maximum and reports them, so a bad Hessian is caught before it is
interpreted as physics.

# Returns
NamedTuple with `infidelity`, `gradient_norm`, `λ_min`, `λ_max`, `psd_violation`
(`-λ_min/λ_max`, zero when PSD) and `ok`.
"""
function hessian_health(fp::FidelityRobustGRAPEProblem, x::Vector{<:Real},
                        H::AbstractMatrix; rtol::Real=1e-6, warn::Bool=true)
    F, g, _, _ = calculate_fidelity_and_derivatives(fp, x)
    λ = eigvals(Symmetric(Matrix(H)))
    λmax = maximum(λ); λmin = minimum(λ)
    viol = λmax > 0 ? max(0.0, -λmin / λmax) : Inf
    ok = viol <= rtol
    if warn && !ok
        @warn "hessian_health: Hessian is not positive semidefinite at a claimed maximum" *
              " (λ_min/λ_max = $(λmin / λmax)). A large ‖dt·H0‖ (e.g. an explicit big" *
              " blockade) breaks the finite-difference step; try the perfect-blockade" *
              " limit or a larger ε_hess." infidelity=1 - F gradient_norm=norm(g)
    end
    return (infidelity=1 - F, gradient_norm=norm(g), λ_min=λmin, λ_max=λmax,
            psd_violation=viol, ok=ok)
end


"""
    apply_control_constraint(H, C, nb) -> (Hc, idx_reduced)

Restrict a Hessian to a linear reparameterization `x = C·y` of the controls, returning
``C^T H C`` and the indices of the additional parameters in the reduced vector.

Used when the physical knobs are fewer than the propagation's parameters -- a composite
gate that plays ONE waveform in each of several segments is the motivating case: the
propagation carries `n_segments × ntimes` steps, but the experiment programs `ntimes`.

`C === nothing` is a no-op. Otherwise `C` must be `N × m` and must leave the additional
parameters alone: its last `nb` columns must be the last `nb` unit vectors, and the
additional-parameter rows of its first `m - nb` columns must vanish. Both are checked --
without them the trailing block of the reduced Hessian is not `x_add` any more, and the
subsequent Schur complement would silently eliminate the wrong thing.
"""
function apply_control_constraint(H::AbstractMatrix,
                                  C::Union{Nothing,AbstractMatrix}, nb::Int)
    N = size(H, 1)
    C === nothing && return (Matrix(H), (N - nb + 1):N)
    size(C, 1) == N || throw(ArgumentError(
        "apply_control_constraint: constraint has $(size(C,1)) rows, Hessian is $N×$N"))
    m = size(C, 2)
    m > nb || throw(ArgumentError("apply_control_constraint: constraint leaves no controls"))
    for j in 1:nb
        col = C[:, m - nb + j]
        e = zeros(eltype(col), N); e[N - nb + j] = 1
        col == e || throw(ArgumentError(
            "apply_control_constraint: column $(m - nb + j) must be unit vector " *
            "$(N - nb + j) so the additional parameters pass through untouched"))
    end
    if nb > 0 && m - nb > 0
        iszero(C[(N - nb + 1):N, 1:(m - nb)]) || throw(ArgumentError(
            "apply_control_constraint: the control columns must not touch the " *
            "additional-parameter rows"))
    end
    return (Matrix(Symmetric(transpose(C) * H * C)), (m - nb + 1):m)
end
