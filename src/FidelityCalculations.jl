using LinearAlgebra, FFTW

"""
    calculate_fidelity_and_derivatives(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real})

Calculate the fidelity between the evolved unitary and target unitary, along with its derivatives.

# Arguments
- `fidelity_problem::FidelityRobustGRAPEProblem`: The fidelity problem definition
- `x::Vector{<:Real}`: The optimization vector containing control parameters and additional parameters

# Returns
A tuple with:
- `F`: The fidelity value
- `F_dx_tot`: Combined derivatives with respect to all control parameters (both main and additional)
- `F_d2err`: Second derivatives with respect to error sources
- `F_d2err_dx_tot`: Combined mixed derivatives for error and all control parameters
"""
function calculate_fidelity_and_derivatives(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real})
    unitary_problem = fidelity_problem.unitary_problem
    ndim = unitary_problem.ndim
    (U, U_dx, U_dx_add, U_derr, U_derr_dx, U_derr_dx_add) = calculate_unitary_and_derivatives(unitary_problem, x)
    ntimes = unitary_problem.ntimes
    nb_additional_param = unitary_problem.nb_additional_param
    nerr = length(unitary_problem.error_sources)
    
    x_main = x[1:end-nb_additional_param]
    nparam = div(length(x_main), ntimes)
    x_main = reshape(x_main, (nparam, ntimes))
    x_add = x[end-nb_additional_param+1:end]

    U0 = fidelity_problem.target_unitary(x_add)
    U0_dx_add = zeros(Complex,ndim,ndim,nb_additional_param)
    x_add_copy = copy(x_add)
    for npa=1:nb_additional_param
        x_add_copy[npa] += unitary_problem.ϵ
        U0_temp = fidelity_problem.target_unitary(x_add_copy)
        U0_dx_add[:,:,npa] = (1/unitary_problem.ϵ) * (U0_temp-U0)
        x_add_copy[npa] = x_add[npa]
    end
    F_dx = zeros(Real, nparam, ntimes)
    F_dx_add = zeros(Real, nb_additional_param)
    F_d2err = zeros(Real, nerr)
    F_d2err_dx = zeros(Real, nparam, ntimes, nerr)
    F_d2err_dx_add = zeros(Real, nb_additional_param, nerr)
    
    P0 = convert.(Complex,fidelity_problem.projector)
    tr_mod(A) = tr(P0*A)
    P = copy(P0)
    P[P .!= 0] .= 1
    D = real(tr(P0))

    
    F = (real(tr_mod(P * U0' * U * P * U' * U0)) + abs(tr_mod(P * U0' * U))^2)/(D*(D+1))

    for nt=1:ntimes
        for np=1:nparam
            F_dx[np,nt] = (
                real(tr_mod(
                    P * U0' * U_dx[:,:,np,nt] * P * U' * U0
                    + P * U0' * U * P * (U_dx[:,:,np,nt])' * U0
                )) + 2*real(conj(tr_mod(P * U0' * U)) * tr_mod(P * U0' * U_dx[:,:,np,nt]))
            )/(D*(D+1))
        end
    end
    
    for npa=1:nb_additional_param
        F_dx_add[npa] = (
            real(tr_mod(
                P * U0' * U_dx_add[:,:,npa] * P * U' * U0
                + P * U0' * U * P * (U_dx_add[:,:,npa])' * U0
                + P * (U0_dx_add[:,:,npa])' * U * P * U' * U0
                + P * U0' * U * P * U' * U0_dx_add[:,:,npa]
            )) + 2*real(conj(tr_mod(P * U0' * U)) * tr_mod(P * U0' * U_dx_add[:,:,npa] + P * (U0_dx_add[:,:,npa])' * U))
        )/(D*(D+1))
    end
    
    for ne=1:nerr
        F_d2err[ne] = 2*(
            real(tr_mod(P * U0' * U_derr[:,:,ne] * P * (U_derr[:,:,ne])' * U0 - P * (U_derr[:,:,ne])' * U_derr[:,:,ne]))
            + abs(tr_mod(P * U0' * U_derr[:,:,ne]))^2
            -D*real(tr_mod(P * (U_derr[:,:,ne])' * U_derr[:,:,ne]))
        )/(D*(D+1))
        
        for nt=1:ntimes
            for np=1:nparam
                F_d2err_dx[np,nt,ne] = 2*(
                    real(tr_mod(
                        P * U0' * U_derr_dx[:,:,np,nt,ne] * P * (U_derr[:,:,ne])' * U0
                        + P * U0' * U_derr[:,:,ne] * P * (U_derr_dx[:,:,np,nt,ne])' * U0
                        - P * (U_derr_dx[:,:,np,nt,ne])' * U_derr[:,:,ne]
                        - P * (U_derr[:,:,ne])' * U_derr_dx[:,:,np,nt,ne]
                    )) + 2*real(conj(tr_mod(P * U0' * U_derr[:,:,ne])) * tr_mod(P * U0' * U_derr_dx[:,:,np,nt,ne]))
                    - D*real(tr_mod(P * (U_derr_dx[:,:,np,nt,ne])' * U_derr[:,:,ne] + P * (U_derr[:,:,ne])' * U_derr_dx[:,:,np,nt,ne]))
                )/(D*(D+1))
            end
        end
        
        for npa=1:nb_additional_param
            F_d2err_dx_add[npa,ne] = 2*(
                real(tr_mod(
                    P * (U0_dx_add[:,:,npa])' * U_derr[:,:,ne] * P * (U_derr[:,:,ne])' * U0
                    + P * U0' * U_derr_dx_add[:,:,npa,ne] * P * (U_derr[:,:,ne])' * U0
                    + P * U0' * U_derr[:,:,ne] * P * (U_derr_dx_add[:,:,npa,ne])' * U0
                    + P * U0' * U_derr[:,:,ne] * P * (U_derr[:,:,ne])' * U0_dx_add[:,:,npa]
                    - P * (U_derr_dx_add[:,:,npa,ne])' * U_derr[:,:,ne]
                    - P * (U_derr[:,:,ne])' * U_derr_dx_add[:,:,npa,ne]
                )) + 2*real(conj(tr_mod(P * U0' * U_derr[:,:,ne])) *
                    tr_mod(P * (U0_dx_add[:,:,npa])' * U_derr[:,:,ne] + P * U0' * U_derr_dx_add[:,:,npa,ne]))
                - D*real(tr_mod(P * (U_derr_dx_add[:,:,npa,ne])' * U_derr[:,:,ne]
                    + P * (U_derr[:,:,ne])' * U_derr_dx_add[:,:,npa,ne]))
            )/(D*(D+1))
        end
    end
    
    F_dx_tot = [reshape(F_dx,nparam*ntimes); F_dx_add]
    F_d2err_dx_tot = cat(reshape(F_d2err_dx,nparam*ntimes,nerr),F_d2err_dx_add;dims=1)
    return (F, F_dx_tot, F_d2err, F_d2err_dx_tot)
end

"""
    optimize_fidelity_and_error_sources(fidelity_problem::FidelityRobustGRAPEProblem, fidelity_parameters::FidelityRobustGRAPEParameters)

Optimizes quantum control pulses while considering robustness against error sources.

This is a high-level wrapper around Optim.optimize that:
1. Combines fidelity optimization with error robustness
2. Applies regularization to the control pulses
3. Handles caching to avoid redundant calculations

# Parameters
- `fidelity_problem::FidelityRobustGRAPEProblem`: Problem definition including the Hamiltonian, target unitary, and error sources
- `fidelity_parameters::FidelityRobustGRAPEParameters`: Optimization parameters including initial values, regularization, and solver configuration

# Returns
- An Optim.OptimizationResults object containing the optimization status, parameters, and diagnostics

# Example
```julia
# Create a problem definition
problem = FidelityRobustGRAPEProblem(...)

# Create optimization parameters
params = FidelityRobustGRAPEParameters(
    x_initial = initial_x,
    regularization_functions = [regularization_cost_phase],
    regularization_coeff1 = [1e-6],
    regularization_coeff2 = [1e-6],
    error_source_coeff = [1e-3],
    iterations = 1000,
    solver_algorithm = LBFGS()
)

# Run the optimization
result = optimize_fidelity_and_error_sources(problem, params)

# Get optimized parameters
optimal_params = Optim.minimizer(result)
```
"""
function optimize_fidelity_and_error_sources(fidelity_problem::FidelityRobustGRAPEProblem, fidelity_parameters::FidelityRobustGRAPEParameters)
    @assert size(fidelity_parameters.error_source_coeff) == size(fidelity_problem.unitary_problem.error_sources)

    nerr = size(fidelity_problem.unitary_problem.error_sources,1)
    ntimes = fidelity_problem.unitary_problem.ntimes
    nb_additional_param = fidelity_problem.unitary_problem.nb_additional_param
    nparam = (size(fidelity_parameters.x_initial,1) - nb_additional_param)÷ntimes

    @assert size(fidelity_parameters.regularization_coeff1,1) == nparam
    @assert size(fidelity_parameters.regularization_coeff2,1) == nparam
    @assert size(fidelity_parameters.regularization_functions,1) == nparam


    function calculate_common!(x, last_x, buffer)
        if x != last_x
            copy!(last_x, x)
            F, F_dx, F_d2err, F_d2err_dx = calculate_fidelity_and_derivatives(fidelity_problem,x)
            buffer[1] = 1-F
            buffer[2:end] = - F_dx
            if size(F_d2err,1) > 0
                buffer[1] += sum(fidelity_parameters.error_source_coeff .* F_d2err .^ 2)
                F_d2err_squared = 2 * sum(reshape(fidelity_parameters.error_source_coeff .* F_d2err,1,nerr) .* F_d2err_dx, dims=2)[:,1]
                buffer[2:end] .+= F_d2err_squared
            end

            reg_costs_grad = zeros(nparam,ntimes)
            reg_costs_tot = zeros(nparam)
            x_main = reshape(x[1:end-nb_additional_param], nparam, ntimes)
            for np=1:nparam
                r1,j1,r2,j2 = fidelity_parameters.regularization_functions[np](x_main[np,:])
                reg_costs_tot[np] = fidelity_parameters.regularization_coeff1[np]*r1+fidelity_parameters.regularization_coeff2[np]*r2
                reg_costs_grad[np,:] = fidelity_parameters.regularization_coeff1[np]*j1+fidelity_parameters.regularization_coeff2[np]*j2
            end
            buffer[1] += sum(reg_costs_tot)
            buffer[2:end-nb_additional_param] += sum(reg_costs_grad,dims=1)[1,:]
        end
    end
    
    function f(x,buffer,last_x)
        calculate_common!(x,last_x,buffer)
        return buffer[1]
    end
    
    function g!(x,stor,buffer,last_x)
        calculate_common!(x, last_x, buffer)
        stor .= buffer[2:end]
    end

    buffer = zeros(typeof(fidelity_parameters.x_initial[1]), size(fidelity_parameters.x_initial,1)+1)
    last_x = similar(fidelity_parameters.x_initial)
    res_optim = optimize(x -> f(x,buffer,last_x), (stor,x) -> g!(x,stor,buffer,last_x), fidelity_parameters.x_initial;
        method = fidelity_parameters.solver_algorithm,
        time_limit = fidelity_parameters.time_limit,
        iterations = fidelity_parameters.iterations,
        fidelity_parameters.additional_parameters...
    )
    return res_optim
end


"""
    calculate_fidelity_response(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real}, normalized_frequencies::Vector{<:Real})

Calculate the frequency-domain fidelity response function for error sources at specified frequencies.

The fidelity response function characterizes how sensitive the quantum gate is to noise at different 
frequencies. This is crucial for understanding robustness against various noise spectra and can guide
the design of control pulses that are specifically robust against the dominant noise frequencies in
a particular experimental setup.

# Parameters
- `fidelity_problem::FidelityRobustGRAPEProblem`: The fidelity problem definition
- `x::Vector{<:Real}`: The control parameters vector
- `normalized_frequencies::Vector{<:Real}`: Normalized frequencies at which to evaluate the response function

# Returns
- A matrix of dimensions (nfreq, nerr) containing the fidelity response function values for each error 
  source at each frequency, where nfreq is the number of frequencies and nerr is the number of error sources.

# Notes
- Frequencies are normalized by Ω (the characteristic frequency of the system)
- The response function is related to the error sensitivity: at zero frequency, the response function
  value is proportional to the static error sensitivity
- Higher values indicate greater sensitivity to noise at that frequency
"""
function calculate_fidelity_response(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real}, normalized_frequencies::Vector{<:Real})
    unitary_problem = fidelity_problem.unitary_problem
    ndim = unitary_problem.ndim
    ntimes = unitary_problem.ntimes
    nerr = size(unitary_problem.error_sources,1)

    nfreq = size(normalized_frequencies,1)
    dt = fidelity_problem.unitary_problem.t0/ntimes
    error_operators_int = calculate_interaction_error_operators(unitary_problem, x)
    P0 = convert.(Complex,fidelity_problem.projector)
    tr_mod(A) = tr(P0*A)
    P = copy(P0)
    P[P .!= 0] .= 1
    D = real(tr(P0))
    time_indices = collect(0:ntimes-1)
    time_indices_reshaped = reshape(time_indices,1,1,ntimes)

    response_fct_ω = zeros(Real,nfreq,nerr)
    for ne=1:nerr
        for nf=1:nfreq
            ω = normalized_frequencies[nf]
            sum_error_freq = sum((@. exp(-im*ω*dt*time_indices_reshaped) * error_operators_int[:,:,:,ne]);dims=3)[:,:,1]
            r = 0.
            for k=1:ntimes
                r += (
                    1. / D * real(exp(im*ω*dt*k) * tr_mod(error_operators_int[:,:,k,ne] * sum_error_freq * P))
                    - 1. / (D*(D+1)) * real(exp(im*ω*dt*k) * tr_mod(error_operators_int[:,:,k,ne] * P * sum_error_freq * P))
                    - 1. / (D*(D+1)) * real(exp(im*ω*dt*k) * tr_mod(error_operators_int[:,:,k,ne] * P) * tr_mod(sum_error_freq * P))
                )
            end
            response_fct_ω[nf,ne] = dt^2 * r
        end
    end
    return response_fct_ω
end

"""
    calculate_fidelity_response_fft(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real}; oversampling::Int = 1)

Calculate the frequency-domain fidelity response function using Fast Fourier Transform (FFT).

This function efficiently computes the frequency-domain fidelity response function for all error sources
using FFT algorithms, which is computationally more efficient than direct calculation at specific
frequencies. The implementation adds zero-padding to achieve frequency oversampling if requested.

# Parameters
- `fidelity_problem::FidelityRobustGRAPEProblem`: The fidelity problem definition
- `x::Vector{<:Real}`: The control parameters vector
- `oversampling::Int = 1`: Oversampling factor to increase the resolution of the frequency response

# Returns
- A tuple containing:
  - `response_fct_ω`: Matrix of dimensions (ntimes*oversampling, nerr) containing the fidelity 
    response function values for each frequency and error source
  - `norm_frequencies`: Vector of normalized frequencies corresponding to the response function values

# Notes
- Frequencies are normalized by the quantum system's characteristic energy scale
- The oversampling parameter allows for higher frequency resolution by zero-padding in the time domain
"""
function calculate_fidelity_response_fft(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real}; oversampling::Int = 1)
    @assert oversampling >= 1
    unitary_problem = fidelity_problem.unitary_problem
    ndim = unitary_problem.ndim
    ntimes = unitary_problem.ntimes
    nerr = size(unitary_problem.error_sources,1)

    dt = fidelity_problem.unitary_problem.t0/ntimes
    error_operators_int = calculate_interaction_error_operators(unitary_problem, x)
    error_operators_int = convert.(ComplexF64,error_operators_int)

    zero_pad = zeros(ComplexF64,ndim,ndim,ntimes*(oversampling-1),nerr)
    error_operators_int = cat(error_operators_int,zero_pad;dims=3)
    # (ndim,ndim,ntimes,nerr)
    P0 = convert.(Complex,fidelity_problem.projector)
    tr_mod(A) = tr(P0*A)
    P = copy(P0)
    P[P .!= 0] .= 1
    D = real(tr(P0))

    response_fct_ω = zeros(Real,ntimes*oversampling,nerr)
    for ne=1:nerr
        p_fft = plan_fft(error_operators_int[:,:,:,ne], (3,))
        error_operators_int_fft = p_fft * error_operators_int[:,:,:,ne]
        p_ifft = plan_ifft(error_operators_int[:,:,:,ne], (3,))
        error_operators_int_ifft = (ntimes*oversampling) * (p_ifft * error_operators_int[:,:,:,ne])

        for nt=1:(ntimes*oversampling)
            response_fct_ω[nt,ne] = dt^2 * (
                1/D * real(tr_mod(error_operators_int_ifft[:,:,nt] * error_operators_int_fft[:,:,nt] * P))
                - 1/(D*(D+1)) * real(tr_mod(error_operators_int_ifft[:,:,nt] * P * error_operators_int_fft[:,:,nt] * P))
                - 1/(D*(D+1)) * real(tr_mod(error_operators_int_ifft[:,:,nt] * P) * tr_mod(error_operators_int_fft[:,:,nt] * P))
            )
        end
    end
    norm_frequencies = (2*π/(ntimes*oversampling*dt)) * collect(0:(ntimes*oversampling-1))
    return response_fct_ω, norm_frequencies
end

"""
    calculate_expectation_values(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real})

Calculate the time-dependent expectation values of error generators during the quantum evolution.

This function computes how the expectation values of the error generators evolve over time
under the optimized control pulse sequence. This provides insight into how errors accumulate
during the quantum gate implementation and can be used to visualize the error sensitivity
profile over time.

# Parameters
- `fidelity_problem::FidelityRobustGRAPEProblem`: The fidelity problem definition
- `x::Vector{<:Real}`: The control parameters vector

# Returns
- A matrix of dimensions (ntimes, nerr) containing the expectation values of each error
  operator at each time step of the evolution.

# Notes
- The calculation uses the cumulative sum of the interaction-picture error operators
- Values are normalized by the dimension of the projector subspace
- The time resolution is determined by the number of time steps in the problem definition
"""
function calculate_expectation_values(fidelity_problem::FidelityRobustGRAPEProblem, x::Vector{<:Real})
    unitary_problem = fidelity_problem.unitary_problem
    ndim = unitary_problem.ndim
    ntimes = unitary_problem.ntimes
    error_operators_int = calculate_interaction_error_operators(unitary_problem, x)
    nerr = size(error_operators_int,4)
    error_operators_cumsum = cumsum(error_operators_int; dims=3)
    exp_values = zeros(Real,ntimes,nerr)
    dt = unitary_problem.t0/ntimes

    P0 = convert.(Complex,fidelity_problem.projector)
    tr_mod(A) = tr(P0*A)
    P = copy(P0)
    P[P .!= 0] .= 1
    D = real(tr(P0))

    for ne=1:nerr
        for nt=1:ntimes
            exp_values[nt,ne] = real(dt*tr_mod(error_operators_cumsum[:,:,nt,ne])/D)
        end
    end
    return exp_values
end