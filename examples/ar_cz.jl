##
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Random
using Optim
using Setfield
using FunctionWrappers: FunctionWrapper

##
Random.seed!(43)
ntimes = 200
t0_to = 7.613
t0_ar = 14.32

# Define Hamiltonian function and wrap it
H0_func(time_step, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
H0 = HamiltonianFunctionWrapper(H0_func)

# Define error Hamiltonians and wrap them
H_amplitude_error_func(time_step, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0(time_step, ϕ, x_add)
H_amplitude_error = ErrorHamiltonianFunctionWrapper(H_amplitude_error_func)

# Define target unitary function and wrap it
cz_func(x_add) = cz_with_1q_phase_symmetric(x_add[1])
cz = UnitaryFunctionWrapper(cz_func)

rydberg_problem = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=t0_to,
        ntimes=ntimes,
        ndim=5,
        H0=H0,
        nb_additional_param=1,
        error_sources=[]
    ),
    collect(Float64,Diagonal([1, 2, 1, 0, 0])),
    cz
)

rydberg_problem_ar = (@set rydberg_problem.unitary_problem.error_sources = [
    ErrorSource(H_amplitude_error)
])
rydberg_problem_ar = (@set rydberg_problem_ar.unitary_problem.t0 = t0_ar)

rydberg_cz_parameters = FidelityRobustGRAPEParameters(
    x_initial = [2*π*0.001*rand(ntimes); 2*π*rand()],
    # Use wrapped regularization function for better performance
    regularization_functions = [regularization_cost_phase_wrapped],
    regularization_coeff1=[1e-7],
    regularization_coeff2=[1e-7],
    error_source_coeff=Vector{Real}(),
    time_limit=40,
    additional_parameters = Dict(
        :show_trace => true,
        :show_every => 10,
        :g_tol => 1e-9
    )
)

rydberg_cz_parameters_ar = (@set rydberg_cz_parameters.error_source_coeff = [1e-4])
rydberg_cz_parameters_ar = (@set rydberg_cz_parameters_ar.time_limit = 300)
rydberg_cz_parameters_ar = (@set rydberg_cz_parameters_ar.regularization_coeff1 = [1e-6])
rydberg_cz_parameters_ar = (@set rydberg_cz_parameters_ar.regularization_coeff2 = [1e-6])

res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, rydberg_cz_parameters)
res_optim_cz_ar = optimize_fidelity_and_error_sources(rydberg_problem_ar, rydberg_cz_parameters_ar)

##

optim_pulse = Optim.minimizer(res_optim_cz)
optim_pulse_ar = Optim.minimizer(res_optim_cz_ar)
using PyPlot
fig, ax = subplots()
ax.plot((0:(ntimes-1)) / ntimes * t0_to, unwrap_phase(optim_pulse[1:ntimes]), label="Time-optimal")
ax.plot((0:(ntimes-1)) / ntimes * t0_ar, unwrap_phase(optim_pulse_ar[1:ntimes]), label="Amplitude-error robust")
ax.set_xlabel("Time (1/Ω)")
ax.set_ylabel("Laser phase (rad)")
ax.legend()
@show fig
##

# Define frequency error Hamiltonian and wrap it
H_frequency_error_func(time_step, ϕ, x_add, δ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, δ) - H0(time_step, ϕ, x_add)
H_frequency_error = ErrorHamiltonianFunctionWrapper(H_frequency_error_func)

rydberg_problem_with_errors = (@set rydberg_problem.unitary_problem.error_sources = [
    ErrorSource(H_amplitude_error),
    ErrorSource(H_frequency_error)
])

rydberg_problem_with_errors_ar = (@set rydberg_problem_ar.unitary_problem.error_sources = [
    ErrorSource(H_amplitude_error),
    ErrorSource(H_frequency_error)
])

F, _, F_d2err, _ = calculate_fidelity_and_derivatives(rydberg_problem_with_errors, optim_pulse)
F_ar, _, F_d2err_ar, _ = calculate_fidelity_and_derivatives(rydberg_problem_with_errors_ar, optim_pulse_ar)

# Define decay operator and wrap it
decay_operator_func(time_step, x, x_add, ϵ) = ϵ*collect(ComplexF64,Diagonal([0, 0, 0, 1, 1]))
decay_operator = ErrorHamiltonianFunctionWrapper(decay_operator_func)

rydberg_problem_with_decay = (@set rydberg_problem.unitary_problem.error_sources = [
    ErrorSource(decay_operator)
])
rydberg_problem_with_decay_ar = (@set rydberg_problem_ar.unitary_problem.error_sources = [
    ErrorSource(decay_operator)
])
rydberg_pop = calculate_expectation_values(rydberg_problem_with_decay, optim_pulse)[end, 1]
rydberg_pop_ar = calculate_expectation_values(rydberg_problem_with_decay_ar, optim_pulse)[end, 1]

println("Time-optimal gate")
println("Sensitivity to amplitude errors: F = 1 - $(-F_d2err[1]/2) × ϵ²")
println("Sensitivity to frequency errors: F = 1 - $(-F_d2err[2]/2) × δ²")
println("Integrated Rydberg population: $(rydberg_pop)/Ω")

println("\nAmplitude-error robust gate")
println("Sensitivity to amplitude errors: F = 1 - $(-F_d2err_ar[1]/2) × ϵ²")
println("Sensitivity to frequency errors: F = 1 - $(-F_d2err_ar[2]/2) × δ²")
println("Integrated Rydberg population: $(rydberg_pop_ar)/Ω")
##

#frequencies = collect(LinRange(0,3,100))
response_fct, frequencies = calculate_fidelity_response_fft(rydberg_problem_with_errors, optim_pulse; oversampling=30)
response_fct_ar, frequencies_ar = calculate_fidelity_response_fft(rydberg_problem_with_errors_ar, optim_pulse_ar; oversampling=30)
@assert isapprox(-F_d2err[1], 2*response_fct[1,1], rtol=1e-3, atol=1e-3)
@assert isapprox(-F_d2err[2], 2*response_fct[1,2], rtol=1e-3, atol=1e-3)
@assert isapprox(-F_d2err_ar[1], 2*response_fct_ar[1,1], rtol=1e-3, atol=1e-3)
@assert isapprox(-F_d2err_ar[2], 2*response_fct_ar[1,2], rtol=1e-3, atol=1e-3)
##
idx_error=1
fig, ax = subplots()
ax.plot(frequencies, 0.25 * response_fct[:,idx_error], label="Time-optimal")
ax.plot(frequencies_ar, 0.25 * response_fct_ar[:,idx_error], label="Amplitude-error robust")
ax.set_xlim(0,5)
ax.set_ylim(0,4)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser intensity noise fidelity response")
ax.set_title("Time-optimal vs amplitude-error robust")
ax.legend()
@show fig

##
idx_error=2
fig, ax = subplots()
ax.plot(frequencies, (2*π)^2 * response_fct[:,idx_error], label="Time-optimal")
ax.plot(frequencies_ar, (2*π)^2 * response_fct_ar[:,idx_error], label="Amplitude-error robust")
ax.set_xlim(0,3.5)
ax.set_ylim(0,500)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser frequency noise fidelity response")
ax.set_title("Time-optimal vs amplitude-error robust")
ax.legend()
@show fig