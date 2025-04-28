using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Test
using Random
using Optim

# Helper function for phase regularization
function regularization_cost_phase(x)
    # Compute phase difference, handling 2π periodicity
    diff_x = diff(x)
    diff_diff_x = diff(diff_x)
    
    # First and second order regularization
    reg1 = sum(sin.(diff_x/2).^2)
    reg2 = sum(sin.(diff_diff_x/2).^2)
    
    # Gradients
    n = length(x)
    jac1 = zeros(n)
    jac2 = zeros(n)
    
    for i in 1:n-1
        if i < n-1
            jac1[i] -= 0.5*sin(diff_x[i])
        end
        if i > 1
            jac1[i] += 0.5*sin(diff_x[i-1])
        end
    end
    
    for i in 1:n
        if i < n-2
            jac2[i] -= 0.5*sin(diff_diff_x[i])
        end
        if i > 1 && i < n-1
            jac2[i] += sin(diff_diff_x[i-1])
        end
        if i > 2
            jac2[i] -= 0.5*sin(diff_diff_x[i-2])
        end
    end
    
    return reg1, jac1, reg2, jac2
end

@testset "RobustGRAPE.jl" begin
    @testset "Error sensitivity gradient validation for CZ gate" begin
        println("[TEST] Validating error sensitivity gradient for CZ gate...")
        
        # Setup parameters exactly matching the test_cz_grad.ipynb notebook
        ntimes = 200
        t0 = 2*π*1.22
        nb_tests = 2  # Number of gradient tests to run
        
        # Define the Hamiltonian, error Hamiltonian, and target gate
        H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])
        Herror(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0(t, ϕ, x_add)
        
        # Create the problem with error source
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0,
                nb_additional_param=1,
                error_sources=[
                    ErrorSource(Herror=Herror)
                ]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Projector onto computational subspace with factor 2 for |01⟩
            target_unitary=cz
        )
        
        # Set random seed for reproducibility
        Random.seed!(42)
        
        # Run multiple tests of the error sensitivity gradient calculation
        for ntest = 1:nb_tests
            # Choose a parameter to test the gradient for
            # For the last test, test the additional parameter, otherwise test a random time step
            if ntest == nb_tests
                idx_test = ntimes + 1
            else
                idx_test = rand(1:ntimes)
            end
            
            # Generate random control parameters
            xs = 2*π*rand(ntimes + 1)
            
            # Use same perturbation size as in the notebook
            ϵ_test = 1e-4
            
            # Calculate error sensitivity and its gradient at initial point
            _, _, F_d2err0, F_d2err0_dx = calculate_fidelity_and_derivatives(rydberg_problem, xs)
            
            # Make a small change to test the gradient
            xs[idx_test] += ϵ_test
            
            # Calculate new error sensitivity
            _, _, F_d2err1, _ = calculate_fidelity_and_derivatives(rydberg_problem, xs)
            
            # Test that numerical approximation (F_d2err1-F_d2err0)/ϵ_test matches the analytical gradient
            @test isapprox(
                (F_d2err1[1] - F_d2err0[1])/ϵ_test, 
                F_d2err0_dx[idx_test, 1], 
                rtol=1e-3, 
                atol=1e-5
            )
        end
    end
    
    @testset "Time-optimal CZ gate" begin
        println("[TEST] Verifying known CZ gate solution fidelity...")
        # This test is based on the time_optimal_cz.ipynb example
        # It verifies that the known solution for a time-optimal CZ gate 
        # achieves high fidelity using RobustGRAPE
        #
        # The parameters used for this pulse were found and described in:
        # Evered, S.J., Bluvstein, D., Kalinowski, M. et al. 
        # High-fidelity parallel entangling gates on a neutral-atom quantum computer.
        # Nature 622, 268–272 (2023)

        # Parameter values for the known solution
        t0 = 2*π*1.22           # Total time duration
        A = 0.7701624           # Amplitude
        ω0 = 0.97525275         # Frequency
        ϕ0 = -0.97449603        # Phase
        δ0 = -0.04319765        # Detuning slope
        ntimes = 1000           # Number of time steps
        θ = 2.0802725844516097  # Single-qubit phase parameter

        # Generate control pulse shape (phase vs time)
        times = LinRange(0, t0, ntimes)
        ϕs = @. (A*cos(ω0*times-ϕ0) + δ0*times)
        xs = [ϕs; θ]  # Complete parameter vector with phase parameter

        # Define the Hamiltonian and target unitary (without amplitude or detuning error)
        H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        # Define the target CZ gate composed with single-qubit phase gates
        # (whose phase is an additional parameter)
        cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])

        # Create the problem definition
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )

        # Calculate fidelity of the known solution
        fidelity = calculate_fidelity_and_derivatives(rydberg_problem, xs)[1]
        
        # Test that fidelity is very high (infidelity < 1e-4)
        @test fidelity > 0.9999
    end
    
    @testset "Pulse optimization and error sensitivity" begin
        println("[TEST] Testing analytical calculation of error sensitivity...")
        # Test pulse optimization and error sensitivity analysis based on time_optimal_cz.ipynb
        
        # Setup problem parameters matching the notebook
        ntimes = 200
        t0 = 2*π*1.22
        
        # Define system and target gate
        H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])
        
        # Create the base problem definition
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )
        
        # Set random seed for reproducibility
        Random.seed!(42)
        
        # Create initial parameters
        initial_ϕ = (2*π*0.001) .* rand(Float64, ntimes)
        initial_x = [initial_ϕ; 2*π*rand()]
        
        # Configure optimization parameters using the FidelityRobustGRAPEParameters
        optimization_params = FidelityRobustGRAPEParameters(
            x_initial = initial_x,
            regularization_functions = [regularization_cost_phase],
            regularization_coeff1 = [1e-6],
            regularization_coeff2 = [1e-6],
            error_source_coeff = Vector{Real}(),  # No error sources in this problem
            iterations = 40,
            solver_algorithm = LBFGS(),
            additional_parameters = Dict(
                :f_abstol => 1e-11,
                :g_tol => 3e-10,
                :show_trace => false
            )
        )
        
        # Run optimization with the wrapper function
        res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, optimization_params)
        
        # Get optimized parameters
        optimized_x = Optim.minimizer(res_optim_cz)
        
        # Calculate baseline fidelity with optimized pulse
        F0, _, _ = calculate_fidelity_and_derivatives(rydberg_problem, optimized_x)
        
        # Test the sensitivity to errors by computing second derivatives
        
        # Define error parameter and Hamiltonians with errors
        error_coeff = rydberg_problem.unitary_problem.ϵ2  # Small error for numerical accuracy
        
        # Error Hamiltonian and perturbed Hamiltonians 
        Herror(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0(t, ϕ, x_add)
        H0_plus_error(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], error_coeff, 0)
        H0_minus_error(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], -error_coeff, 0)
        
        # Create problems with fixed positive/negative errors
        rydberg_problem_plus_error = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0_plus_error,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )
        
        rydberg_problem_minus_error = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0_minus_error,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )
        
        # Create problem with error source for analytical derivative calculation
        robust_problem_with_error = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0,
                nb_additional_param=1,
                error_sources=[
                    ErrorSource(Herror=Herror)
                ]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )
        
        # Calculate fidelities with errors
        F0_plus_error, _, _ = calculate_fidelity_and_derivatives(rydberg_problem_plus_error, optimized_x)
        F0_minus_error, _, _ = calculate_fidelity_and_derivatives(rydberg_problem_minus_error, optimized_x)
        
        # Calculate fidelity and derivatives with error source
        F1, _, F1_d2err = calculate_fidelity_and_derivatives(robust_problem_with_error, optimized_x)
        
        # Compute numerical second derivative using finite difference
        numerical_d2F = (F0_plus_error + F0_minus_error - 2*F0) / error_coeff^2
        
        # Compare numerical and analytical derivatives using the same criteria as in the notebook
        @test isapprox(numerical_d2F, F1_d2err[1], rtol=1e-3, atol=1e-2)
    end
    
    @testset "Fidelity gradient validation" begin
        println("[TEST] Validating gradient calculation accuracy...")
        # This test validates that the analytical gradient calculation in 
        # calculate_fidelity_and_derivatives is accurate by comparing with 
        # numerical finite difference approximation
        
        # Test setup
        t0 = 2*π*1.22                # Arbitrary time duration
        ntimes = 50             # Number of time steps (reduced for test efficiency)
        nb_tests = 5            # Number of gradient tests to run
        
        # Define system model
        H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])
        
        # Create the problem definition
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )
        
        # Set random seed for reproducibility
        Random.seed!(42)
        
        # Run multiple tests of gradient calculation
        for ntest=1:nb_tests
            # Choose a parameter to test the gradient for
            # For the last test, use the additional parameter
            if ntest < nb_tests
                idx_test = rand(1:ntimes)
            else
                idx_test = ntimes+1
            end
            
            # Generate random control parameters
            xs = 2*π*rand(ntimes+1)
            
            # Calculate fidelity and its gradient at the initial point
            F0, grad_F0, _ = calculate_fidelity_and_derivatives(rydberg_problem, xs)
            
            # Make a small change to test the gradient
            xs[idx_test] += rydberg_problem.unitary_problem.ϵ
            
            # Calculate new fidelity
            F1, _, _ = calculate_fidelity_and_derivatives(rydberg_problem, xs)
            
            # Test that numerical approximation (F1-F0)/ϵ matches the analytical gradient
            @test isapprox(
                (F1-F0)/rydberg_problem.unitary_problem.ϵ, 
                grad_F0[idx_test], 
                rtol=1e-3, 
                atol=1e-3
            )
        end
    end
    
    @testset "Gradient-based pulse optimization" begin
        println("[TEST] Testing gradient-based optimization for pulse finding...")
        # This test implements gradient-based optimization using the wrapper function
        # to find a high-fidelity control pulse for a CZ gate
        
        # Setup problem with reduced size for test efficiency
        ntimes = 200
        t0 = 2*π*1.22
        
        # Define system and target gate
        H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])
        
        # Create the problem definition
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz
        )
        
        # Set random seed for reproducibility
        Random.seed!(42)
        
        # Create initial parameters (small random pulse + random phase)
        initial_ϕ = (2*π*0.001) .* rand(Float64, ntimes)
        initial_x = [initial_ϕ; 2*π*rand()]
        
        # Create the optimization parameters structure
        optimization_params = FidelityRobustGRAPEParameters(
            x_initial = initial_x,
            regularization_functions = [regularization_cost_phase],
            regularization_coeff1 = [1e-6],
            regularization_coeff2 = [1e-6],
            error_source_coeff = Vector{Real}(),  # No error sources in this problem
            iterations = 40,
            solver_algorithm = LBFGS(),
            additional_parameters = Dict(
                :f_abstol => 1e-11,
                :g_tol => 3e-10,
                :show_trace => false
            )
        )
        
        # Run optimization with the wrapper function
        res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, optimization_params)

        # Get optimized parameters and calculate final fidelity
        optimized_x = Optim.minimizer(res_optim_cz)
        final_fidelity = calculate_fidelity_and_derivatives(rydberg_problem, optimized_x)[1]
        final_infidelity = 1.0 - final_fidelity
        
        # Verify optimization achieved high fidelity
        @test final_infidelity < 1e-6
    end
    
    @testset "Reduced vs Full Hamiltonian Error Sensitivity" begin
        println("[TEST] Comparing error sensitivity between reduced and full Hamiltonians...")
        
        # Use the time-optimal pulse parameters from previous test
        ntimes = 500  # Match the notebook
        t0 = 7.613    # Match the notebook
        
        # Re-load the optimal pulse for the test
        # Setup problem with reduced size for test efficiency
        H0_sym(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        cz_sym(x_add) = cz_with_1q_phase_symmetric(x_add[1])
        
        # Create the problem definition
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0_sym,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz_sym
        )
        
        # Set random seed for reproducibility
        Random.seed!(43)  # Match the notebook
        
        # Create initial parameters for optimization
        initial_ϕ = (2*π*0.001) .* rand(Float64, ntimes)
        initial_x = [initial_ϕ; 2*π*rand()]
        
        # Create the optimization parameters structure
        optimization_params = FidelityRobustGRAPEParameters(
            x_initial = initial_x,
            regularization_functions = [regularization_cost_phase],
            regularization_coeff1 = [1e-7],  # Match the notebook
            regularization_coeff2 = [1e-7],  # Match the notebook
            error_source_coeff = Vector{Real}(),
            iterations = 40,
            solver_algorithm = LBFGS(),
            additional_parameters = Dict(
                :f_abstol => 1e-11,
                :g_tol => 1e-9,  # Match the notebook
                :show_trace => false
            )
        )
        
        # Run optimization to get the optimized pulse
        res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, optimization_params)
        optim_pulse = Optim.minimizer(res_optim_cz)
        
        # Now check the sensitivity to errors for both Hamiltonian models
        
        # Define the Hamiltonians for the full (non-symmetric) model
        H1_full(t, ϕ, x_add) = rydberg_hamiltonian_full_blockaded(ϕ[1], 0, 0)
        H_intensity_error_full(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_full_blockaded(ϕ[1], ϵ, 0) - H1_full(t, ϕ, x_add)
        H_frequency_error_full(t, ϕ, x_add, δ) = rydberg_hamiltonian_full_blockaded(ϕ[1], 0, δ) - H1_full(t, ϕ, x_add)
        cz_full(x_add) = cz_with_1q_phase_full(x_add[1]; rydberg_dimension=3)
        
        # Create the problem with full Hamiltonian and both error sources
        rydberg_problem_full = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=7,
                H0=H1_full,
                nb_additional_param=1,
                error_sources=[
                    ErrorSource(Herror=H_intensity_error_full),
                    ErrorSource(Herror=H_frequency_error_full)
                ]
            ),
            projector=Diagonal([1,1,1,1,0,0,0]),  # Projector for the full model
            target_unitary=cz_full
        )
        
        # Define the error Hamiltonians for the symmetric model
        H_intensity_error_sym(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0_sym(t, ϕ, x_add)
        H_frequency_error_sym(t, ϕ, x_add, δ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, δ) - H0_sym(t, ϕ, x_add)
        
        # Create the problem with symmetric Hamiltonian and both error sources
        rydberg_problem_sym = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0_sym,
                nb_additional_param=1,
                error_sources=[
                    ErrorSource(Herror=H_intensity_error_sym),
                    ErrorSource(Herror=H_frequency_error_sym)
                ]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz_sym
        )
        
        # Calculate fidelity and error sensitivities for both models
        F_full, _, F_d2err_full, _ = calculate_fidelity_and_derivatives(rydberg_problem_full, optim_pulse)
        F_sym, _, F_d2err_sym, _ = calculate_fidelity_and_derivatives(rydberg_problem_sym, optim_pulse)
        
        # Test that the error sensitivities match between the two models
        # This verifies that the symmetric model correctly captures the error behavior
        @test isapprox(F_d2err_sym[1], F_d2err_full[1], rtol=1e-3, atol=1e-3)
        @test isapprox(F_d2err_sym[2], F_d2err_full[2], rtol=1e-3, atol=1e-3)
        
        # Print the sensitivities for informational purposes
        println("Intensity error sensitivity: $(F_d2err_sym[1])")
        println("Frequency error sensitivity: $(F_d2err_sym[2])")
    end
    
    @testset "Fidelity Response vs Error Sensitivity" begin
        println("[TEST] Verifying fidelity response matches error sensitivity...")
        
        # Use the time-optimal pulse from previous tests
        ntimes = 500  # Match the notebook
        t0 = 7.613    # Match the notebook
        
        # Define the Hamiltonians and setup the problem
        H0_sym(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
        cz_sym(x_add) = cz_with_1q_phase_symmetric(x_add[1])
        
        # Set random seed for reproducibility
        Random.seed!(43)  # Match the notebook
        
        # Create initial parameters for the optimization
        initial_ϕ = (2*π*0.001) .* rand(Float64, ntimes)
        initial_x = [initial_ϕ; 2*π*rand()]
        
        # Create the problem definition for optimization
        rydberg_problem = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0_sym,
                nb_additional_param=1,
                error_sources=[]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz_sym
        )
        
        # Optimization parameters
        optimization_params = FidelityRobustGRAPEParameters(
            x_initial = initial_x,
            regularization_functions = [regularization_cost_phase],
            regularization_coeff1 = [1e-7],  # Match the notebook
            regularization_coeff2 = [1e-7],  # Match the notebook
            error_source_coeff = Vector{Real}(),
            iterations = 40,
            solver_algorithm = LBFGS(),
            additional_parameters = Dict(
                :f_abstol => 1e-11,
                :g_tol => 1e-9,  # Match the notebook
                :show_trace => false
            )
        )
        
        # Run optimization to get the optimized pulse
        res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, optimization_params)
        optim_pulse = Optim.minimizer(res_optim_cz)
        
        # Define error Hamiltonians for the error sensitivity calculation
        H_intensity_error_sym(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0_sym(t, ϕ, x_add)
        H_frequency_error_sym(t, ϕ, x_add, δ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, δ) - H0_sym(t, ϕ, x_add)
        
        # Create the problem with symmetric Hamiltonian and both error sources
        rydberg_problem_with_errors = FidelityRobustGRAPEProblem(
            unitary_problem=UnitaryRobustGRAPEProblem(
                t0=t0,
                ntimes=ntimes,
                ndim=5,
                H0=H0_sym,
                nb_additional_param=1,
                error_sources=[
                    ErrorSource(Herror=H_intensity_error_sym),
                    ErrorSource(Herror=H_frequency_error_sym)
                ]
            ),
            projector=Diagonal([1,2,1,0,0]),  # Corrected with factor 2 for |01⟩
            target_unitary=cz_sym
        )
        
        # Calculate fidelity and its derivatives with respect to errors
        _, _, F_d2err, _ = calculate_fidelity_and_derivatives(rydberg_problem_with_errors, optim_pulse)
        
        # Calculate the fidelity response function at different frequencies
        frequencies = collect(LinRange(0, 3, 100))
        response_fct = calculate_fidelity_response(rydberg_problem_with_errors, optim_pulse, frequencies)
        
        # Compare the error sensitivity to the first term of the fidelity response function
        # The error sensitivity should match the first term (frequency = 0) of the response function
        @test isapprox(-F_d2err[1], 2*response_fct[1,1], rtol=1e-3, atol=1e-3)
        @test isapprox(-F_d2err[2], 2*response_fct[1,2], rtol=1e-3, atol=1e-3)
        
        # Print the results for informational purposes
        println("Intensity error sensitivity: $(F_d2err[1]), Response function: $(response_fct[1,1])")
        println("Frequency error sensitivity: $(F_d2err[2]), Response function: $(response_fct[1,2])")
    end
end