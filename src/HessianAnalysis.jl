using LinearAlgebra

"""
    calculate_fidelity_hessian(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real}; ε_hess::Real=1e-5, verbose::Bool=false)

Compute the Hessian matrix ``-\\partial^2 F / \\partial x_i \\partial x_j`` analytically.

The main-parameter block is computed analytically by extending the piecewise unitary derivative
infrastructure. For parameters at different time steps ``t_i \\neq t_j``, the second derivative
of the unitary ``d^2U/(dx_i dx_j)`` is a product of the interaction-picture first derivatives
(no extra finite differences needed). For same-time-step pairs ``t_i = t_j``, a central finite
difference of the matrix exponential provides ``d^2 U_t / (dx_i dx_j)``.

Additional parameters (e.g. target gate phases) are handled by finite-differencing the
analytical gradient from `calculate_fidelity_and_derivatives`.

The sign convention is chosen so that at a fidelity maximum (F ≈ 1), the returned matrix is
positive semi-definite:

```math
F(x + \\delta x) \\approx 1 - \\frac{1}{2} \\delta x^T H \\delta x
```

# Arguments
- `fidelity_problem`: The fidelity problem definition
- `x`: Parameter vector (should be near a fidelity maximum)
- `ε_hess`: Step size for finite differences of the matrix exponential (default: 1e-5)
- `verbose`: Print progress (default: false)

# Returns
- Symmetric matrix of size `(length(x), length(x))`
"""
function calculate_fidelity_hessian(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real};
                                     ε_hess::Real=1e-5, verbose::Bool=false)
    unitary_problem = fidelity_problem.unitary_problem
    ntimes = unitary_problem.ntimes
    ndim = unitary_problem.ndim
    nb_add = unitary_problem.nb_additional_param

    x_main = x[1:end-nb_add]
    @assert mod(length(x_main), ntimes) == 0
    nparam = div(length(x_main), ntimes)
    N_main = nparam * ntimes
    N_total = N_main + nb_add
    x_main_2d = reshape(copy(x_main), (nparam, ntimes))
    x_add = x[end-nb_add+1:end]

    dt = unitary_problem.t0 / ntimes

    # ================================================================
    # Forward pass: compute interaction-picture first and second derivatives
    # ================================================================
    cum_evo = Matrix{ComplexF64}(I, ndim, ndim)
    infimU_dx = zeros(ComplexF64, ndim, ndim, nparam, ntimes)
    infimU_d2x = zeros(ComplexF64, ndim, ndim, nparam, nparam, ntimes)

    x_copy = zeros(nparam)

    for nt in 1:ntimes
        if verbose && (nt % 50 == 0 || nt == 1)
            println("  Forward pass: step $nt / $ntimes")
        end

        H0_base = unitary_problem.H0(nt, @view(x_main_2d[:, nt]), x_add)
        infim_evo = exp(-im * dt * H0_base)
        old_cum_evo = copy(cum_evo)
        cum_evo = infim_evo * cum_evo
        cum_evo_inv = cum_evo'

        infim_p = zeros(ComplexF64, ndim, ndim, nparam)
        infim_m = zeros(ComplexF64, ndim, ndim, nparam)

        x_copy .= @view(x_main_2d[:, nt])
        for np in 1:nparam
            x_copy[np] = x_main_2d[np, nt] + ε_hess
            infim_p[:, :, np] = exp(-im * dt * unitary_problem.H0(nt, x_copy, x_add))
            x_copy[np] = x_main_2d[np, nt] - ε_hess
            infim_m[:, :, np] = exp(-im * dt * unitary_problem.H0(nt, x_copy, x_add))
            x_copy[np] = x_main_2d[np, nt]

            # Central FD first derivative
            infimU_dx[:, :, np, nt] = cum_evo_inv * (1 / (2 * ε_hess)) *
                (infim_p[:, :, np] - infim_m[:, :, np]) * old_cum_evo

            # Central FD second derivative (same parameter)
            infimU_d2x[:, :, np, np, nt] = cum_evo_inv * (1 / ε_hess^2) *
                (infim_p[:, :, np] + infim_m[:, :, np] - 2 * infim_evo) * old_cum_evo
        end

        # Cross second derivatives at the same time step (different parameters)
        for np1 in 1:nparam, np2 in (np1+1):nparam
            x_copy .= @view(x_main_2d[:, nt])
            x_copy[np1] += ε_hess; x_copy[np2] += ε_hess
            infim_pp = exp(-im * dt * unitary_problem.H0(nt, x_copy, x_add))

            x_copy .= @view(x_main_2d[:, nt])
            x_copy[np1] += ε_hess; x_copy[np2] -= ε_hess
            infim_pm = exp(-im * dt * unitary_problem.H0(nt, x_copy, x_add))

            x_copy .= @view(x_main_2d[:, nt])
            x_copy[np1] -= ε_hess; x_copy[np2] += ε_hess
            infim_mp = exp(-im * dt * unitary_problem.H0(nt, x_copy, x_add))

            x_copy .= @view(x_main_2d[:, nt])
            x_copy[np1] -= ε_hess; x_copy[np2] -= ε_hess
            infim_mm = exp(-im * dt * unitary_problem.H0(nt, x_copy, x_add))

            infimU_d2x[:, :, np1, np2, nt] = cum_evo_inv * (1 / (4 * ε_hess^2)) *
                (infim_pp - infim_pm - infim_mp + infim_mm) * old_cum_evo
            infimU_d2x[:, :, np2, np1, nt] = infimU_d2x[:, :, np1, np2, nt]
        end
    end

    U_final = cum_evo

    # ================================================================
    # Fidelity quantities
    # ================================================================
    U0 = fidelity_problem.target_unitary(x_add)
    P0 = convert.(ComplexF64, fidelity_problem.projector)
    P = copy(P0)
    P[P .!= 0] .= 1
    D = real(tr(P0))

    W = U0' * U_final
    T = tr(P0 * W)
    PW_dag = P * W'

    # Precompute W_i = U₀†·U_i and T_i = tr(P₀·W_i) for each main parameter
    W_i = zeros(ComplexF64, ndim, ndim, N_main)
    T_i = zeros(ComplexF64, N_main)
    for nt in 1:ntimes, np in 1:nparam
        idx = (nt - 1) * nparam + np
        W_i[:, :, idx] = W * infimU_dx[:, :, np, nt]
        T_i[idx] = tr(P0 * @view(W_i[:, :, idx]))
    end

    # ================================================================
    # Assemble main-main Hessian block
    #
    # d²F/(dx_i dx_j) = [2Re(tr(P₀ W_ij P W†)) + 2Re(tr(P₀ W_i P W_j†))
    #                   + 2Re(conj(T_j)·T_i) + 2Re(conj(T)·T_ij)] / (D(D+1))
    #
    # where W_ij = U₀†·d²U/(dx_i dx_j) and:
    #   t_i > t_j: W_ij = W · Ã_i · Ã_j       (W = U₀†U, Ã = interaction-picture deriv)
    #   t_j > t_i: W_ij = W · Ã_j · Ã_i
    #   t_i = t_j: W_ij = W · infimU_d2x[:,:,np_i,np_j,t]
    # ================================================================
    if verbose
        println("  Assembling Hessian ($N_main × $N_main main block)...")
    end

    H_main = zeros(N_main, N_main)

    for idx_i in 1:N_main
        np_i = ((idx_i - 1) % nparam) + 1
        nt_i = ((idx_i - 1) ÷ nparam) + 1
        Wi = @view(W_i[:, :, idx_i])
        Ti = T_i[idx_i]

        for idx_j in idx_i:N_main
            np_j = ((idx_j - 1) % nparam) + 1
            nt_j = ((idx_j - 1) ÷ nparam) + 1
            Wj = @view(W_i[:, :, idx_j])
            Tj = T_i[idx_j]

            # Terms that only need W_i, W_j
            term1 = 2 * real(tr(P0 * Wi * P * Wj')) + 2 * real(conj(Tj) * Ti)

            # Terms that need W_ij
            if nt_i > nt_j
                W_ij = W * infimU_dx[:, :, np_i, nt_i] * infimU_dx[:, :, np_j, nt_j]
            elseif nt_j > nt_i
                W_ij = W * infimU_dx[:, :, np_j, nt_j] * infimU_dx[:, :, np_i, nt_i]
            else
                W_ij = W * infimU_d2x[:, :, np_i, np_j, nt_i]
            end

            T_ij = tr(P0 * W_ij)
            term2 = 2 * real(tr(P0 * W_ij * PW_dag)) + 2 * real(conj(T) * T_ij)

            H_main[idx_i, idx_j] = (term1 + term2) / (D * (D + 1))
            if idx_i != idx_j
                H_main[idx_j, idx_i] = H_main[idx_i, idx_j]
            end
        end
    end

    # ================================================================
    # Additional parameters via FD of the analytical gradient
    # ================================================================
    H_full = zeros(N_total, N_total)
    H_full[1:N_main, 1:N_main] = H_main

    if nb_add > 0
        _, grad0, _, _ = calculate_fidelity_and_derivatives(fidelity_problem, x)

        for npa in 1:nb_add
            if verbose
                println("  FD for additional parameter $npa / $nb_add")
            end
            x_pert = copy(x)
            x_pert[N_main + npa] += ε_hess
            _, grad_pert, _, _ = calculate_fidelity_and_derivatives(fidelity_problem, x_pert)

            col = (grad_pert - grad0) / ε_hess
            H_full[:, N_main + npa] = col
            H_full[N_main + npa, :] = col
        end

        # Symmetrize the add-add block
        for npa1 in 1:nb_add, npa2 in (npa1+1):nb_add
            avg = (H_full[N_main+npa1, N_main+npa2] + H_full[N_main+npa2, N_main+npa1]) / 2
            H_full[N_main+npa1, N_main+npa2] = avg
            H_full[N_main+npa2, N_main+npa1] = avg
        end
    end

    return -H_full
end


"""
    calculate_principal_parameters(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real};
                                    n_eigenvalues::Int=10, ε_hess::Real=1e-5, verbose::Bool=false)

Compute the principal parameter functions of a fidelity landscape.

At a fidelity maximum, the Hessian ``-\\partial^2 F / \\partial x_i \\partial x_j`` characterizes
how the fidelity degrades under parameter perturbations. Its eigenvectors define "principal
parameter functions" — the directions in parameter space that the fidelity is most sensitive to.

For problems with many control parameters but few physical constraints (e.g., Rydberg gates
with phase-only control), only a handful of eigenvalues will be significantly nonzero.
The fidelity near the optimum expands as:

```math
F(x + \\delta x) \\approx 1 - \\sum_k \\frac{\\lambda_k}{2} (v_k \\cdot \\delta x)^2
```

# Arguments
- `fidelity_problem`: The fidelity problem definition
- `x`: Optimized parameter vector (should satisfy F ≈ 1)
- `n_eigenvalues`: Maximum number of eigenvalue/eigenvector pairs to return (default: 10)
- `ε_hess`: Step size for Hessian computation (default: 1e-5)
- `verbose`: Print progress (default: false)

# Returns
A named tuple with fields:
- `eigenvalues`: Vector of the largest eigenvalues (sorted descending), length ≤ `n_eigenvalues`
- `eigenvectors`: Matrix whose columns are the corresponding eigenvectors
- `hessian`: The full Hessian matrix ``-\\partial^2 F / \\partial x_i \\partial x_j``
"""
function calculate_principal_parameters(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real};
                                         n_eigenvalues::Int=10, ε_hess::Real=1e-5, verbose::Bool=false)
    H = calculate_fidelity_hessian(fidelity_problem, x; ε_hess=ε_hess, verbose=verbose)

    eig = eigen(Symmetric(H))
    λ = eig.values
    V = eig.vectors

    perm = sortperm(λ; rev=true)
    λ = λ[perm]
    V = V[:, perm]

    n_keep = min(n_eigenvalues, length(λ))

    return (eigenvalues=λ[1:n_keep], eigenvectors=V[:, 1:n_keep], hessian=H)
end
