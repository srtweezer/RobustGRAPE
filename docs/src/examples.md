# Examples

## Time-Optimal CZ Gate

This example shows how to optimize a CZ gate in a Rydberg atom system with minimal gate time.

```julia
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Optim
using Random

# Set random seed for reproducibility
Random.seed!(43)

# Define parameters
ntimes = 500  # Number of time steps
t0 = 7.613    # Total evolution time (in units of 1/Ω)

# Define Hamiltonian and target operation
H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])

# Create the optimization problem
rydberg_problem = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=t0,
        ntimes=ntimes,
        ndim=5,
        H0=H0,
        nb_additional_param=1,
        error_sources=[]
    ),
    Diagonal([1, 2, 1, 0, 0]),  # Target state projection
    cz                          # Target operation
)

# Configure optimization parameters
rydberg_cz_parameters = FidelityRobustGRAPEParameters(
    x_initial = [2*π*0.001*rand(ntimes); 2*π*rand()],  # Initial pulse + phase parameter
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

# Run optimization
res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, rydberg_cz_parameters)

# Extract optimized pulse
optim_pulse = Optim.minimizer(res_optim_cz)
```

## Robust Control Pulse Design

This example demonstrates designing control pulses that are robust against specific noise sources like intensity and frequency errors.

```julia
# Define error Hamiltonians
H_intensity_error(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0(t, ϕ, x_add)
H_frequency_error(t, ϕ, x_add, δ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, δ) - H0(t, ϕ, x_add)

# Create a GRAPE problem with error sources
rydberg_problem_with_errors = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=t0,
        ntimes=ntimes,
        ndim=5,
        H0=H0,
        nb_additional_param=1,
        error_sources=[
            ErrorSource(H_intensity_error),
            ErrorSource(H_frequency_error)
        ]
    ),
    Diagonal([1, 2, 1, 0, 0]),
    cz
)

# Analyze sensitivity to errors with previously optimized pulse
F, _, F_d2err, _ = calculate_fidelity_and_derivatives(rydberg_problem_with_errors, optim_pulse)

println("Infidelity: $(1-F)")
println("Sensitivity to intensity errors: F = 1 - $(-F_d2err[1]/2) × ϵ²")
println("Sensitivity to frequency errors: F = 1 - $(-F_d2err[2]/2) × δ²")
```

## Frequency Domain Analysis

This example shows how to analyze the frequency response of optimized control pulses to understand their robustness characteristics.

```julia
# Calculate frequency response with oversampling for better resolution
response_fct, frequencies = calculate_fidelity_response_fft(
    rydberg_problem_with_errors, 
    optim_pulse; 
    oversampling=30
)

# Verify that DC component matches our direct calculation
@assert isapprox(-F_d2err[1], 2*response_fct[1,1], rtol=1e-3, atol=1e-3)
@assert isapprox(-F_d2err[2], 2*response_fct[1,2], rtol=1e-3, atol=1e-3)

# Plot the response to intensity noise (example plotting with PyPlot)
using PyPlot
fig, ax = subplots()
ax.plot(frequencies, 0.25 * response_fct[:,1])
ax.set_xlim(0, 3)
ax.set_ylim(0, 1.5)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser intensity noise fidelity response")
ax.set_title("Time-optimal gate")

# Plot the response to frequency noise
fig, ax = subplots()
ax.plot(frequencies, (2*π)^2 * response_fct[:,2])
ax.set_xlim(0, 3)
ax.set_ylim(0, 250)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser frequency noise fidelity response")
ax.set_title("Time-optimal gate")
```

## Additional Analysis: Rydberg Population

You can also analyze other properties of the optimized pulse:

```julia
# Calculate Rydberg state population during gate execution
rydberg_problem_with_decay = deepcopy(rydberg_problem_with_errors)
decay_operator(t, x, x_add, ϵ) = ϵ*collect(Diagonal([0, 0, 0, 1, 1]))
rydberg_problem_with_decay = (@set rydberg_problem_with_decay.unitary_problem.error_sources = [
    ErrorSource(decay_operator)
])
rydberg_pop = calculate_expectation_values(rydberg_problem_with_decay, optim_pulse)[end, 1]
println("Integrated Rydberg population: $(rydberg_pop)/Ω")
```

For complete examples, see the `examples/` directory in the package repository.