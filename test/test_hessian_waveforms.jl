#
# Tests for HessianWaveforms.jl, plus the two-quadrature coverage that
# `Fidelity Hessian validation` in runtests.jl does not provide.
#
# THE GAP THIS CLOSES. The existing Hessian test runs a problem with nparam = 1, so the
# same-time-step CROSS-quadrature block of `calculate_fidelity_hessian` -- the 4-point
# stencil over (np1, np2) with np1 != np2 -- has never been exercised. That block is
# exactly the phase-to-modulus coupling, i.e. everything a polar-coordinate principal
# waveform is made of.
#
# WHY MOST CHECKS SIT AT AN OPTIMUM. `F_d2err` is assembled from the FIRST-order U_derr
# alone (see FidelityCalculations.jl); the d2U/derr2 terms are absent because they cancel
# against unitarity only when U reproduces the target. It is therefore the true d2F/derr2
# only where F = 1. Measured on a random pulse of this very problem it is off by 27.7x
# (amplitude) and 6.0x (frequency), and at the optimum it agrees to 5 digits. Any test
# that compares a Hessian contraction against F_d2err must be run at F = 1 or it is
# measuring the wrong thing -- confidently, and with no error raised.
#
# The identities used below:
#
#   gauge        shifting every phase by a constant conjugates H by exp(i c n_r), which is
#                trivial on the computational block, so F is EXACTLY invariant (verified
#                to 1e-16) and the Hessian must annihilate that direction. Holds at ANY
#                pulse, which makes it the sharpest available probe of the Hessian itself.
#   detuning     a phase ramp is a constant detuning in a rotating frame, so the two
#                parameterizations are the same function of the ramp slope and all their
#                derivatives agree. Limited to O(dt^2) only by sampling the ramp at step
#                midpoints.
#   amplitude    a uniform relative-modulus offset IS the intensity error, exactly.
#
# The known-bad case for the conditioning guard is built so the physics is provably
# untouched: H0 -> H0 + c*I multiplies U by a global phase, which cancels out of the
# fidelity exactly, while making ||dt*H0|| enormous. Any change is then pure numerics.
#

@testset "Hessian waveforms" begin

    # ── shared two-quadrature setup: p = [phase, relative modulus] ──────────────────
    # `rydberg_hamiltonian_symmetric_blockaded(ϕ, ϵ, δ)` already carries the (1+ϵ)
    # amplitude scaling, so the second control is the modulus deviation with no extra
    # algebra, and the perfect blockade keeps ||dt*H0|| ~ 1.
    NDIM = 5
    PROJ = collect(Diagonal([1, 2, 1, 0, 0]))
    EH = 1e-4          # ε_hess: 1e-4 leaves the gauge residual at 8e-9 here, 1e-5 at 8e-7
    H0_2q(nt, p, xa) = rydberg_hamiltonian_symmetric_blockaded(p[1], p[2], 0.0)
    Hfrq(nt, p, xa, δ) = rydberg_hamiltonian_symmetric_blockaded(p[1], p[2], δ) -
                         rydberg_hamiltonian_symmetric_blockaded(p[1], p[2], 0.0)
    Hamp(nt, p, xa, ε) = rydberg_hamiltonian_symmetric_blockaded(p[1], p[2] + ε, 0.0) -
                         rydberg_hamiltonian_symmetric_blockaded(p[1], p[2], 0.0)

    function make_problem(t0, ntimes; sources=ErrorSource[], shift=0.0)
        H0 = shift == 0.0 ? H0_2q :
             (nt, p, xa) -> H0_2q(nt, p, xa) + shift * Matrix(I, NDIM, NDIM)
        FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0, ntimes=ntimes, ndim=NDIM, H0=H0,
                nb_additional_param=1, error_sources=sources),
            projector=PROJ,
            target_unitary=xa -> cz_with_1q_phase_symmetric(xa[1]))
    end

    pack(ϕ, a, θ) = vcat(vec(permutedims(hcat(ϕ, a))), θ)
    iphase(nt) = 1:2:2nt
    imod(nt) = 2:2:2nt

    # ── a genuine optimum, built once and shared ───────────────────────────────────
    # Optimized with PHASE ONLY and then embedded at zero modulus deviation: that is the
    # physical situation (fixed envelope, phase-modulated drive), and because F = 1 is the
    # global maximum over ALL parameters the point is still a maximum of the two-quadrature
    # problem, so its Hessian is positive semidefinite and the modulus quadrature has zero
    # gradient despite never having been optimized.
    T_OPT = 2π * 1.22
    NT_OPT = 50
    xopt, fp_opt, fpe_opt = let ntimes = NT_OPT, t0 = T_OPT
        fp1 = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0, ntimes=ntimes, ndim=NDIM,
                H0=(nt, p, xa) -> rydberg_hamiltonian_symmetric_blockaded(p[1], 0.0, 0.0),
                nb_additional_param=1, error_sources=ErrorSource[]),
            projector=PROJ, target_unitary=xa -> cz_with_1q_phase_symmetric(xa[1]))
        Random.seed!(42)
        pars = FidelityRobustGRAPEParameters(
            x_initial=[(2π * 0.001) .* rand(ntimes); 2π * rand()],
            regularization_functions=[regularization_cost_phase],
            regularization_coeff1=[0.0], regularization_coeff2=[0.0],
            error_source_coeff=Vector{Real}(), iterations=3000,
            solver_algorithm=LBFGS(), additional_parameters=Dict(:show_trace => false))
        z = Optim.minimizer(optimize_fidelity_and_error_sources(fp1, pars))
        (pack(z[1:ntimes], zeros(ntimes), z[end]),
         make_problem(t0, ntimes),
         make_problem(t0, ntimes; sources=[ErrorSource(Hamp), ErrorSource(Hfrq)]))
    end
    println("  shared optimum: 1-F = $(1 - calculate_fidelity_and_derivatives(fp_opt, xopt)[1])")
    @test calculate_fidelity_and_derivatives(fp_opt, xopt)[1] > 1 - 1e-12
    # the modulus quadrature was never optimized, yet its gradient must vanish
    @test norm(calculate_fidelity_and_derivatives(fp_opt, xopt)[2][imod(NT_OPT)]) < 1e-7

    # ── 1. two-quadrature Hessian against numerical second derivatives ──────────────
    @testset "nparam=2 Hessian vs finite differences" begin
        ntimes = 14
        fp = make_problem(2π * 1.22, ntimes)
        Random.seed!(2026)
        x = pack(2π * rand(ntimes), 0.15 * randn(ntimes), 2π * rand())
        N = length(x)

        H = calculate_fidelity_hessian(fp, x; ε_hess=EH)
        @test isapprox(H, H', atol=1e-10)

        F(v) = calculate_fidelity_and_derivatives(fp, v)[1]
        function d2num(i, j; h=1e-4)
            xpp = copy(x); xpp[i] += h; xpp[j] += h
            xpm = copy(x); xpm[i] += h; xpm[j] -= h
            xmp = copy(x); xmp[i] -= h; xmp[j] += h
            xmm = copy(x); xmm[i] -= h; xmm[j] -= h
            -(F(xpp) - F(xpm) - F(xmp) + F(xmm)) / (4h^2)
        end

        # The block the existing suite cannot reach: SAME time step, DIFFERENT quadrature.
        worst = 0.0
        for nt in (1, 5, ntimes)
            i = (nt - 1) * 2 + 1        # phase at step nt
            j = (nt - 1) * 2 + 2        # modulus at the same step
            @test isapprox(H[i, j], d2num(i, j); rtol=1e-3, atol=1e-7)
            @test isapprox(H[j, j], d2num(j, j); rtol=1e-3, atol=1e-7)
            worst = max(worst, abs(H[i, j] - d2num(i, j)))
        end
        println("  same-step cross-quadrature block: worst |analytic - FD| = $worst")

        Random.seed!(99)
        for _ in 1:6
            i, j = rand(1:N), rand(1:N)
            @test isapprox(H[i, j], d2num(i, j); rtol=1e-2, atol=1e-6)
        end
    end

    # ── 2. exact identities, at the optimum ────────────────────────────────────────
    @testset "gauge, detuning and amplitude identities" begin
        ntimes = NT_OPT
        _, _, F_d2err, _ = calculate_fidelity_and_derivatives(fpe_opt, xopt)
        H = calculate_fidelity_hessian(fp_opt, xopt; ε_hess=EH)

        # (a) uniform phase offset is an exact symmetry -> exact null direction
        vgauge = zeros(length(xopt)); vgauge[iphase(ntimes)] .= 1.0
        println("  gauge   ||H·1_phase||  = $(norm(H * vgauge))")
        @test norm(H * vgauge) < 1e-6

        # (b) uniform relative modulus IS the amplitude error, with no discretization
        vamp = zeros(length(xopt)); vamp[imod(ntimes)] .= 1.0
        println("  amp     ΣH_aa = $(dot(vamp, H, vamp))   -d²F/dε² = $(-F_d2err[1])")
        @test isapprox(dot(vamp, H, vamp), -F_d2err[1]; rtol=1e-5)

        # (c) a phase ramp is a constant detuning in a rotating frame; sampling it at step
        #     midpoints is what limits this to O(dt²)
        dt = T_OPT / ntimes
        vramp = zeros(length(xopt))
        vramp[iphase(ntimes)] .= [(k - 0.5) * dt for k in 1:ntimes]
        println("  detune  tᵀHt  = $(dot(vramp, H, vramp))   -d²F/dδ² = $(-F_d2err[2])")
        @test isapprox(dot(vramp, H, vramp), -F_d2err[2]; rtol=1e-2)

        # (d) at a maximum the Hessian is positive semidefinite and LOW RANK
        λ = eigvals(Symmetric(H))
        println("  spectrum: λ_max = $(maximum(λ)), λ_min = $(minimum(λ)), " *
                "rank(>1e-4·λ₁) = $(count(λ .> 1e-4 * maximum(λ))) of $(length(λ))")
        @test minimum(λ) > -1e-6 * maximum(λ)
        @test count(λ .> 1e-4 * maximum(λ)) < length(λ) ÷ 4
    end

    # ── 3. conditioning guard, verified against known-good AND known-bad ────────────
    @testset "conditioning guard" begin
        ntimes = 20
        t0 = 2π * 1.22
        Random.seed!(11)
        x = pack(2π * rand(ntimes), zeros(ntimes), 2π * rand())
        vgauge = zeros(length(x)); vgauge[iphase(ntimes)] .= 1.0

        fp_ok = make_problem(t0, ntimes)
        fp_bad = make_problem(t0, ntimes; shift=1e5)

        # H0 + c·I is a global phase on U, so the fidelity is analytically unchanged; what
        # little it moves is itself the roundoff this test is about.
        dF = abs(calculate_fidelity_and_derivatives(fp_ok, x)[1] -
                 calculate_fidelity_and_derivatives(fp_bad, x)[1])
        println("  |ΔF| from the c·I shift = $dF  (analytically zero)")
        @test dF < 1e-8

        g_ok = norm(calculate_fidelity_hessian(fp_ok, x; ε_hess=EH) * vgauge)
        g_bad = norm(calculate_fidelity_hessian(fp_bad, x; ε_hess=EH) * vgauge)
        println("  ||H·1_phase||: well conditioned $g_ok , ‖dt·H0‖~4e4 $g_bad " *
                "(ratio $(g_bad / g_ok))")
        @test g_ok < 1e-6                # known good
        @test g_bad / g_ok > 1e3         # known bad -- the detector must actually fire
    end

    # ── 4. Schur complement == curvature of the θ-maximized fidelity ────────────────
    @testset "schur_eliminate vs max over additional parameters" begin
        ntimes = 12
        fp = make_problem(2π * 1.22, ntimes)
        Random.seed!(31)
        u0 = vec(permutedims(hcat(2π * rand(ntimes), 0.1 * randn(ntimes))))
        N = 2ntimes + 1

        Fθ(u, θ) = calculate_fidelity_and_derivatives(fp, vcat(u, θ))[1]
        θgrid = range(0, 2π, length=721)
        θ0 = θgrid[argmax([Fθ(u0, θ) for θ in θgrid])]
        brent(u) = Optim.optimize(θ -> -Fθ(u, θ), θ0 - 0.4, θ0 + 0.4, Optim.Brent();
                                  abs_tol=1e-14)
        G(u) = -Optim.minimum(brent(u))                 # max_θ F, computed independently
        θstar = Optim.minimizer(brent(u0))

        H = calculate_fidelity_hessian(fp, vcat(u0, θstar); ε_hess=EH)
        Hs = schur_eliminate(H, [N])
        @test size(Hs) == (2ntimes, 2ntimes)
        @test isapprox(Hs, Hs', atol=1e-10)

        function d2G(i, j; h=2e-3)
            up = copy(u0); up[i] += h; up[j] += h
            pm = copy(u0); pm[i] += h; pm[j] -= h
            mp = copy(u0); mp[i] -= h; mp[j] += h
            mm = copy(u0); mm[i] -= h; mm[j] -= h
            -(G(up) - G(pm) - G(mp) + G(mm)) / (4h^2)
        end
        worst = 0.0
        for (i, j) in ((1, 1), (2, 2), (1, 2), (3, 8), (2ntimes - 1, 2ntimes))
            @test isapprox(Hs[i, j], d2G(i, j); rtol=2e-2, atol=1e-5)
            worst = max(worst, abs(Hs[i, j] - d2G(i, j)))
        end
        println("  Schur vs FD of max_θ F: worst |diff| = $worst")

        # guards against a "Schur complement" that silently returns H[keep, keep]
        @test !isapprox(Hs, H[1:2ntimes, 1:2ntimes]; rtol=1e-3)
        @test schur_eliminate(H, Int[]) == Matrix(H)
        @test_throws ErrorException schur_eliminate(zeros(3, 3), [3])
    end

    # ── 5. optimize_additional_parameters ──────────────────────────────────────────
    @testset "optimize_additional_parameters" begin
        ntimes = 12
        fp = make_problem(2π * 1.22, ntimes)
        Random.seed!(5)
        u0 = vec(permutedims(hcat(2π * rand(ntimes), 0.1 * randn(ntimes))))
        Fθ(θ) = calculate_fidelity_and_derivatives(fp, vcat(u0, θ))[1]
        θgrid = range(0, 2π, length=721)
        θ0 = θgrid[argmax([Fθ(θ) for θ in θgrid])]

        xo = optimize_additional_parameters(fp, vcat(u0, θ0))
        g = calculate_fidelity_and_derivatives(fp, xo)[2][end]
        # The floor is set by RobustGRAPE's own forward-difference of the TARGET unitary
        # w.r.t. x_add (UnitaryRobustGRAPEProblem.ϵ = 1e-8), so dF/dθ cannot be resolved
        # below ~1e-9 no matter how many Newton steps are taken.
        println("  x_add Newton: |dF/dθ| = $(abs(g)), Δθ = $(xo[end] - θ0)")
        @test abs(g) < 1e-8
        @test Fθ(xo[end]) >= Fθ(θ0) - 1e-12

        fp0 = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=1.0, ntimes=4, ndim=NDIM, H0=H0_2q,
                nb_additional_param=0, error_sources=ErrorSource[]),
            projector=PROJ, target_unitary=xa -> cz_with_1q_phase_symmetric(0.0))
        v = collect(1.0:8.0)
        @test optimize_additional_parameters(fp0, v) == v
    end

    # ── 6. shift_error_source is the same physics as the error source ───────────────
    @testset "shift_error_source" begin
        F0 = calculate_fidelity_and_derivatives(fpe_opt, xopt)[1]
        _, _, F_d2err, _ = calculate_fidelity_and_derivatives(fpe_opt, xopt)
        for (k, name) in ((1, "amplitude"), (2, "frequency"))
            h = 1e-3
            Fp = calculate_fidelity_and_derivatives(shift_error_source(fpe_opt, k, +h), xopt)[1]
            Fm = calculate_fidelity_and_derivatives(shift_error_source(fpe_opt, k, -h), xopt)[1]
            num = (Fp - 2F0 + Fm) / h^2
            println("  $name: FD of shifted problem $num vs analytic $(F_d2err[k])")
            @test isapprox(num, F_d2err[k]; rtol=1e-4)
        end
        @test_throws ArgumentError shift_error_source(fpe_opt, 5, 0.1)
        # the shifted problem must not silently keep the error sources
        @test isempty(shift_error_source(fpe_opt, 1, 0.1).unitary_problem.error_sources)
        @test length(shift_error_source(fpe_opt, 1, 0.1;
                                        keep_sources=true).unitary_problem.error_sources) == 2
    end

    # ── 7. robustness Hessian: σ must cancel out ───────────────────────────────────
    @testset "calculate_robustness_hessian" begin
        ntimes = NT_OPT
        fpe = make_problem(T_OPT, ntimes; sources=[ErrorSource(Hfrq)])
        r1 = calculate_robustness_hessian(fpe, xopt; σ=0.08, ε_hess=EH)
        r2 = calculate_robustness_hessian(fpe, xopt; σ=0.04, ε_hess=EH)
        @test size(r1.K) == (2ntimes, 2ntimes)
        @test isapprox(r1.K, r1.K', atol=1e-8)

        rel = norm(r1.K - r2.K) / norm(r1.K)
        println("  ||K(σ=.08) - K(σ=.04)|| / ||K|| = $rel")
        @test rel < 2e-2
        # the eigenVECTORS are what the figure shows, so check the leading subspace too
        e1 = eigen(Symmetric(r1.K)); e2 = eigen(Symmetric(r2.K))
        ov = abs(dot(e1.vectors[:, end], e2.vectors[:, end]))
        println("  leading eigenvector overlap = $ov ; λ₁ = $(e1.values[end])")
        @test ov > 0.999
        # NOT tested here: positive semidefiniteness. K = 4G'G only where the gate is
        # robust; the shared optimum is the TIME-OPTIMAL gate, whose detuning sensitivity
        # is 5.7, so K legitimately has negative eigenvalues. Instead, validate K against
        # a direct fourth-order finite difference of the same quantity, which is what it
        # claims to be regardless of robustness.
        let σ = 0.05, hu = 1e-2
            Fmax(u, s) = begin                # max_θ F at static detuning s
                fps = s == 0 ? fpe : shift_error_source(fpe, 1, s)
                Fθ(θ) = calculate_fidelity_and_derivatives(fps, vcat(u, θ))[1]
                θ0 = xopt[end]
                -Optim.minimum(Optim.optimize(θ -> -Fθ(θ), θ0 - 0.5, θ0 + 0.5,
                                              Optim.Brent(); abs_tol=1e-14))
            end
            # A(u) = -d²F̃/dε², by the SAME ±σ average the implementation uses, so with the
            # same σ the two agree exactly rather than only to O(σ²).
            A(u) = -(Fmax(u, σ) + Fmax(u, -σ) - 2 * Fmax(u, 0)) / σ^2
            u0 = xopt[1:end-1]
            function d2A(i, j)
                pp = copy(u0); pp[i] += hu; pp[j] += hu
                pm = copy(u0); pm[i] += hu; pm[j] -= hu
                mp = copy(u0); mp[i] -= hu; mp[j] += hu
                mm = copy(u0); mm[i] -= hu; mm[j] -= hu
                (A(pp) - A(pm) - A(mp) + A(mm)) / (4hu^2)
            end
            r = calculate_robustness_hessian(fpe, xopt; σ=σ, ε_hess=EH)
            worst = 0.0
            for (i, j) in ((1, 1), (2, 2), (1, 2), (7, 7), (5, 20), (2ntimes - 1, 2ntimes))
                num = d2A(i, j)
                @test isapprox(r.K[i, j], num; rtol=3e-2, atol=2e-3)
                worst = max(worst, abs(r.K[i, j] - num))
            end
            println("  K vs direct 4th-order FD: worst |diff| = $worst " *
                    "(‖K‖_∞ = $(maximum(abs, r.K)))")
        end

        # σ = 0 block must agree with the plain Schur-eliminated Hessian at the SAME x_add
        xr = vcat(xopt[1:end-1], r1.x_add[1])
        H = calculate_fidelity_hessian(fpe, xr; ε_hess=EH)
        @test isapprox(r1.fidelity_hessian, schur_eliminate(H, [length(xr)]); rtol=1e-6)
        @test isapprox(r1.robustness_cost_hessian,
                       r1.fidelity_hessian + r1.sigma^2 / 2 * r1.K; atol=1e-8)
    end

    # ── 8. principal_waveforms ─────────────────────────────────────────────────────
    @testset "principal_waveforms" begin
        ntimes = NT_OPT
        H = schur_eliminate(calculate_fidelity_hessian(fp_opt, xopt; ε_hess=EH),
                            [length(xopt)])
        pw = principal_waveforms(H, ntimes, 2; n_max=6)
        @test length(pw.waveforms) == 6
        @test all(size(w) == (2, ntimes) for w in pw.waveforms)
        @test issorted(pw.eigenvalues; rev=true)
        @test isapprox(pw.eigenvectors' * pw.eigenvectors, I(6), atol=1e-10)
        # the reshape must match the package's (nparam, ntimes) layout
        @test pw.waveforms[1][1, 3] == pw.eigenvectors[(3 - 1) * 2 + 1, 1]
        @test pw.waveforms[1][2, 3] == pw.eigenvectors[(3 - 1) * 2 + 2, 1]
        # sign convention is deterministic, not whatever LAPACK returned
        @test pw.eigenvectors == principal_waveforms(H, ntimes, 2; n_max=6).eigenvectors
        for k in 1:6
            v = pw.eigenvectors[:, k]
            @test v[argmax(abs.(v))] > 0
        end
        @test_throws ArgumentError principal_waveforms(H, ntimes, 3)
        println("  rank at rtol=1e-4: $(pw.rank) of $(2ntimes);  λ = " *
                "$(round.(pw.eigenvalues; sigdigits=3))")
        @test pw.rank <= 12
    end

    # ── 9. hessian_health, at a real optimum ───────────────────────────────────────
    @testset "hessian_health" begin
        h_ok = hessian_health(fp_opt, xopt,
                              calculate_fidelity_hessian(fp_opt, xopt; ε_hess=EH))
        println("  at the optimum: 1-F = $(h_ok.infidelity), |∇F| = $(h_ok.gradient_norm), " *
                "psd_violation = $(h_ok.psd_violation)")
        @test h_ok.ok
        @test h_ok.gradient_norm < 1e-6

        # same physics, ruined conditioning
        fp_bad = make_problem(T_OPT, NT_OPT; shift=1e5)
        h_bad = hessian_health(fp_bad, xopt,
                               calculate_fidelity_hessian(fp_bad, xopt; ε_hess=EH);
                               warn=false)
        println("  with ‖dt·H0‖~1e5: psd_violation = $(h_bad.psd_violation)")
        @test !h_bad.ok
        @test h_bad.psd_violation > 1e3 * max(h_ok.psd_violation, 1e-12)
    end
end
