##
using Pkg
Pkg.resolve()
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Optim
using Random
using Setfield

##
Random.seed!(43)
ntimes = 500
t0 = 7.613
H0(time_step,ϕ,x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1],0,0)
cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])

rydberg_problem = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=t0,
        ntimes=ntimes,
        ndim=5,
        H0=H0,
        nb_additional_param=1,
        error_sources=[]
    ),
    collect(Diagonal([1,2,1,0,0])),
    cz
)

rydberg_cz_parameters = FidelityRobustGRAPEParameters(
    x_initial = [2*π*0.001*rand(ntimes); 2*π*rand()],
    regularization_functions = [regularization_cost_phase],
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

res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem,rydberg_cz_parameters)


##

optim_pulse = Optim.minimizer(res_optim_cz)
using PyPlot
fig, ax = subplots()
ax.plot((0:(ntimes-1)) / ntimes * t0, unwrap_phase(optim_pulse[1:ntimes]))
ax.set_title("Time-optimal Rydberg pulse")
ax.set_xlabel("Time (1/Ω)")
ax.set_ylabel("Laser phase (rad)")
@show fig
##

H_amplitude_error(time_step,ϕ,x_add,ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1],ϵ,0) - H0(time_step,ϕ,x_add)
H_frequency_error(time_step,ϕ,x_add,δ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1],0,δ) - H0(time_step,ϕ,x_add)


rydberg_problem_with_errors = (@set rydberg_problem.unitary_problem.error_sources = [
    ErrorSource(H_amplitude_error),
    ErrorSource(H_frequency_error)
])

F, _, F_d2err, _ = calculate_fidelity_and_derivatives(rydberg_problem_with_errors,optim_pulse)
decay_operator(time_step,x,x_add,ϵ) = ϵ*collect(Diagonal([0,0,0,1,1]))
rydberg_problem_with_decay = (@set rydberg_problem.unitary_problem.error_sources = [
    ErrorSource(decay_operator)
])
rydberg_pop = calculate_expectation_values(rydberg_problem_with_decay, optim_pulse)[end,1]
println("Infidelity: $(1-F)")
println("Sensitivity to amplitude errors: F = 1 - $(-F_d2err[1]/2) × ϵ²")
println("Sensitivity to frequency errors: F = 1 - $(-F_d2err[2]/2) × δ²")
println("Integrated Rydberg population: $(rydberg_pop)/Ω")
##

#frequencies = collect(LinRange(0,3,100))
response_fct, frequencies = calculate_fidelity_response_fft(rydberg_problem_with_errors,optim_pulse; oversampling=30)
@assert isapprox(-F_d2err[1], 2*response_fct[1,1],rtol=1e-3,atol=1e-3)
@assert isapprox(-F_d2err[2], 2*response_fct[1,2],rtol=1e-3,atol=1e-3)
##
idx_error=1
fig, ax = subplots()
ax.plot(frequencies,0.25 * response_fct[:,idx_error])
ax.set_xlim(0,3)
ax.set_ylim(0,1.5)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser intensity noise fidelity response")
ax.set_title("Time-optimal gate")
@show fig

##
idx_error=2
fig, ax = subplots()
ax.plot(frequencies,(2*π)^2 * response_fct[:,idx_error] )
ax.set_xlim(0,3)
ax.set_ylim(0,250)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser frequency noise fidelity response")
ax.set_title("Time-optimal gate")
@show fig

