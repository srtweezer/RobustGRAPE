# Examples

## Time-Optimal CZ Gate

This example demonstrates how to optimize a CZ gate in a Rydberg atom system with minimal gate time. The optimization is based on the GRAPE algorithm which divides the control pulse into discrete time steps and optimizes the control parameters at each step.

### Problem Setup and Optimization

First, we import the necessary packages and set up the problem:

```julia
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Optim
using Random
using Setfield

# Set random seed for reproducibility
Random.seed!(43)

# Define optimization parameters
ntimes = 500  # Number of time steps
t0 = 7.613    # Total evolution time (in units of 1/Ω)

# Define Hamiltonian and target operation
# Hamiltonian function takes time t, control parameter ϕ and additional parameters
H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)

# Target operation (CZ gate) with an additional phase parameter
cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])
```

Here, `rydberg_hamiltonian_symmetric_blockaded` represents the Hamiltonian of a Rydberg atom system in a blockade configuration. The function takes the laser phase `ϕ[1]` as its first argument, with additional parameters for intensity and detuning errors (set to 0 here).

Next, we create the optimization problem:

```julia
# Create the optimization problem
rydberg_problem = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=t0,             # Total evolution time
        ntimes=ntimes,     # Number of time steps
        ndim=5,            # Dimension of the Hilbert space (5 for the Rydberg system)
        H0=H0,             # Hamiltonian function
        nb_additional_param=1,  # One additional parameter (phase)
        error_sources=[]   # No error sources for basic optimization
    ),
    Diagonal([1, 2, 1, 0, 0]),  # Target state projection (computational subspace)
    cz                          # Target operation
)
```

The problem setup involves:
- A `UnitaryRobustGRAPEProblem` that defines the physical system and dynamics
- A `FidelityRobustGRAPEProblem` that adds the target operation and subspace projection
- The projection `Diagonal([1, 2, 1, 0, 0])` restricts the optimization to the computational subspace

Now we configure the optimization parameters:

```julia
# Configure optimization parameters
rydberg_cz_parameters = FidelityRobustGRAPEParameters(
    # Initial control pulse: small random values for the time steps + random phase
    x_initial = [2*π*0.001*rand(ntimes); 2*π*rand()],
    
    # Use phase regularization to ensure smooth pulses
    regularization_functions = [regularization_cost_phase],
    regularization_coeff1=[1e-7],  # First derivative regularization weight
    regularization_coeff2=[1e-7],  # Second derivative regularization weight
    
    error_source_coeff=Vector{Real}(),  # No error sources yet
    time_limit=40,  # Maximum optimization time in seconds
    
    # Additional parameters for the optimizer
    additional_parameters = Dict(
        :show_trace => true,    # Show optimization progress
        :show_every => 10,      # Display every 10 iterations
        :g_tol => 1e-9          # Gradient tolerance for convergence
    )
)
```

The `regularization_cost_phase` function is particularly important for phase-based control, as it ensures smooth transitions in the phase values, which is crucial for experimental implementation.

Finally, we run the optimization:

```julia
# Run optimization
res_optim_cz = optimize_fidelity_and_error_sources(rydberg_problem, rydberg_cz_parameters)

# Extract optimized pulse
optim_pulse = Optim.minimizer(res_optim_cz)

# Evaluate the final fidelity
final_fidelity, _, _, _ = calculate_fidelity_and_derivatives(rydberg_problem, optim_pulse)
println("Final fidelity: $(final_fidelity)")
```

### Visualizing the Optimized Pulse

We can visualize the optimized pulse using PyPlot:

```julia
using PyPlot

# Plot the phase profile
fig, ax = subplots()
# Use unwrap_phase to ensure smooth phase representation
ax.plot((1:ntimes) / ntimes * t0, unwrap_phase(optim_pulse[1:ntimes]))
ax.set_title("Time-optimal Rydberg pulse")
ax.set_xlabel("Time (1/Ω)")
ax.set_ylabel("Laser phase (rad)")
```

## Robust Control Pulse Design

A key feature of RobustGRAPE is the ability to analyze and optimize control pulses for robustness against various sources of noise and experimental imperfections. This example demonstrates how to evaluate the sensitivity of an optimized quantum gate to common error sources.

### Defining Error Hamiltonians

First, we need to define how different error sources affect our system Hamiltonian:

```julia
# Define error Hamiltonians as deviations from the ideal Hamiltonian

# Intensity error: variation in the Rabi frequency (laser power)
H_intensity_error(t, ϕ, x_add, ϵ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], ϵ, 0) - H0(t, ϕ, x_add)

# Frequency error: variation in the laser detuning
H_frequency_error(t, ϕ, x_add, δ) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, δ) - H0(t, ϕ, x_add)
```

Each error Hamiltonian represents how the system dynamics change when a specific parameter deviates from its nominal value:
- The intensity error (ϵ) represents fluctuations in the laser power or coupling strength
- The frequency error (δ) represents fluctuations in the laser frequency or detuning

### Creating a Problem with Error Sources

Next, we create a problem definition that includes these error sources:

```julia
# Create a GRAPE problem with error sources
rydberg_problem_with_errors = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=t0,
        ntimes=ntimes,
        ndim=5,
        H0=H0,
        nb_additional_param=1,
        error_sources=[
            ErrorSource(H_intensity_error),  # First error source: intensity fluctuations
            ErrorSource(H_frequency_error)   # Second error source: frequency fluctuations
        ]
    ),
    Diagonal([1, 2, 1, 0, 0]),  # Target state projection
    cz                          # Target operation
)
```

### Analyzing Error Sensitivity

Now we can analyze how sensitive our previously optimized pulse is to these error sources:

```julia
# Analyze sensitivity to errors with previously optimized pulse
F, _, F_d2err, _ = calculate_fidelity_and_derivatives(rydberg_problem_with_errors, optim_pulse)

# The F_d2err vector contains the second derivatives of fidelity with respect to error parameters
# These values tell us how quickly the fidelity decreases with small errors

println("Infidelity: $(1-F)")
println("Sensitivity to intensity errors: F = 1 - $(-F_d2err[1]/2) × ϵ²")
println("Sensitivity to frequency errors: F = 1 - $(-F_d2err[2]/2) × δ²")
```

For small errors, the infidelity (1-F) scales quadratically with the error amplitude. The coefficients `-F_d2err[i]/2` indicate how sensitive the gate is to each error source - smaller values mean greater robustness.

### Designing Robust Pulses

To optimize for robustness, we can include error sensitivity in the optimization objective:

```julia
# Create parameters for robust optimization
robust_parameters = FidelityRobustGRAPEParameters(
    x_initial = optim_pulse,  # Start from previously optimized pulse
    regularization_functions = [regularization_cost_phase],
    regularization_coeff1=[1e-7],
    regularization_coeff2=[1e-7],
    error_source_coeff=[0.5, 0.5],  # Non-zero weights for error sensitivities
    time_limit=60,
    additional_parameters = Dict(
        :show_trace => true,
        :show_every => 10,
        :g_tol => 1e-9
    )
)

# Run robust optimization
# res_robust = optimize_fidelity_and_error_sources(rydberg_problem_with_errors, robust_parameters)
```

Setting non-zero `error_source_coeff` values instructs the optimizer to minimize both the gate infidelity and the sensitivity to errors, creating a control pulse that maintains high fidelity even in the presence of experimental imperfections.

## Frequency Domain Analysis

RobustGRAPE provides powerful tools for analyzing the frequency-dependent response of quantum gates to noise. This frequency domain analysis is crucial for understanding how a gate responds to noise at different timescales, which can inform both gate design and experimental implementation.

### Computing the Frequency Response

The frequency response calculation transforms the time-dependent sensitivity of a quantum gate into the frequency domain:

```julia
# Calculate frequency response with oversampling for better resolution
response_fct, frequencies = calculate_fidelity_response_fft(
    rydberg_problem_with_errors, 
    optim_pulse; 
    oversampling=30  # Higher values provide smoother frequency resolution
)

# Verify that DC component matches our direct calculation
# The DC component (zero frequency) corresponds to static parameter shifts
@assert isapprox(-F_d2err[1], 2*response_fct[1,1], rtol=1e-3, atol=1e-3)
@assert isapprox(-F_d2err[2], 2*response_fct[1,2], rtol=1e-3, atol=1e-3)
```

The `calculate_fidelity_response_fft` function computes how the gate fidelity responds to time-dependent noise at different frequencies. The returned matrices contain:
- `response_fct`: A matrix where each column corresponds to an error source, and each row to a frequency
- `frequencies`: The corresponding frequency values (in units of Ω/2π)

### Visualizing the Intensity Noise Response

Let's plot the response to intensity noise:

```julia
using PyPlot
fig, ax = subplots()
ax.plot(frequencies, 0.25 * response_fct[:,1])
ax.set_xlim(0, 3)
ax.set_ylim(0, 1.5)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser intensity noise fidelity response")
ax.set_title("Time-optimal gate response to intensity noise")
```

This plot shows how sensitive the gate is to intensity fluctuations at different frequencies. Key insights:
- Higher values indicate greater sensitivity to noise at that frequency
- The DC component (f=0) corresponds to static offsets in laser intensity
- The frequency axis is normalized to the Rabi frequency (Ω)
- Peaks in the response indicate frequencies where the gate is particularly vulnerable to noise

### Visualizing the Frequency Noise Response

Similarly, we can plot the response to laser frequency noise:

```julia
fig, ax = subplots()
ax.plot(frequencies, (2*π)^2 * response_fct[:,2])
ax.set_xlim(0, 3)
ax.set_ylim(0, 250)
ax.set_xlabel("2πf/Ω")
ax.set_ylabel("Laser frequency noise fidelity response")
ax.set_title("Time-optimal gate response to frequency noise")
```

The scaling factor `(2*π)^2` adjusts the response to represent the sensitivity to frequency fluctuations in standard units, where the detuning is measured in radians/second.

### Applications of Frequency Analysis

Frequency response analysis has several important applications:
1. **Filter Design**: Design filters to suppress noise at frequencies where the gate is most sensitive
2. **Pulse Optimization**: Modify control pulses to reduce sensitivity at problematic frequencies
3. **Error Budgeting**: Allocate error budget across different noise sources based on their impact
4. **Hardware Requirements**: Determine spectral noise requirements for experimental equipment

## Additional Analysis: Rydberg Population

Beyond fidelity and error sensitivity, RobustGRAPE allows you to analyze other important properties of quantum control pulses. For Rydberg atom systems, one crucial metric is the total time spent in the Rydberg state, as these states are susceptible to decoherence and decay.

### Measuring Rydberg State Population

We can calculate the integrated Rydberg state population during gate execution:

```julia
# Create a specialized problem for measuring Rydberg population
rydberg_problem_with_decay = deepcopy(rydberg_problem_with_errors)

# Define an operator that detects population in Rydberg states
# The [0,0,0,1,1] pattern targets the Rydberg states in our 5-level system
decay_operator(t, x, x_add, ϵ) = ϵ*collect(Diagonal([0, 0, 0, 1, 1]))

# Update the problem to use this detection operator
rydberg_problem_with_decay = (@set rydberg_problem_with_decay.unitary_problem.error_sources = [
    ErrorSource(decay_operator)
])

# Calculate expectation values throughout the gate evolution
# The result is a matrix of values for each time step and each error source
rydberg_pops = calculate_expectation_values(rydberg_problem_with_decay, optim_pulse)

# Extract the final integrated Rydberg population (first error source, end of evolution)
rydberg_pop = rydberg_pops[end, 1]
println("Integrated Rydberg population: $(rydberg_pop)/Ω")
```

The integrated Rydberg population is a dimensionless quantity that, when multiplied by the decay rate, gives the total probability of decay during the gate operation. Lower values indicate a more resilient gate against decoherence.

### Visualizing State Evolution

You can also visualize how different states evolve during the gate operation:

```julia
# Calculate state evolution at multiple time points
time_points = LinRange(0, t0, 100)
state_evolution = calculate_state_evolution(rydberg_problem, optim_pulse, time_points)

# Plot Rydberg state population vs time for an initial state
initial_state = [1, 0, 0, 0, 0]  # Example: starting in |00⟩ state
populations = [abs2(state_evolution[i,4,1]) + abs2(state_evolution[i,5,1]) for i in 1:length(time_points)]

# Plot the result
using PyPlot
fig, ax = subplots()
ax.plot(time_points, populations)
ax.set_xlabel("Time (1/Ω)")
ax.set_ylabel("Rydberg state population")
ax.set_title("Rydberg population during gate execution")
```

By analyzing and optimizing these additional properties, you can design quantum gates that not only achieve high fidelity but are also robust against practical experimental constraints like decoherence and decay.

For complete working examples, see the `examples/` directory in the package repository.