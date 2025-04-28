# Getting Started

## Installation

To install RobustGRAPE.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("RobustGRAPE")  # From the Julia registry once published
# OR
Pkg.add(url="https://github.com/USERNAME/RobustGRAPE.jl")  # Directly from GitHub
```

## Basic Usage

Here's a simple example of using RobustGRAPE.jl to optimize a quantum gate:

```julia
using RobustGRAPE
using LinearAlgebra

# Define system parameters
ntimes = 100              # Number of time steps
total_time = 1.0         # Total evolution time
ndim = 4                 # Dimension of Hilbert space

# Define Hamiltonian function
function H0(t, controls, additional_params)
    # Example: A simple two-qubit system with controls
    H = zeros(ComplexF64, ndim, ndim)
    
    # Add free evolution terms
    # ...
    
    # Add control terms
    # ...
    
    return H
end

# Define error sources
error_sources = [
    ErrorSource(function(t, controls, add_params, eps)
        # Example error Hamiltonian
        return eps * randn() * I(ndim)
    end)
]

# Create problem definition
unitary_problem = UnitaryRobustGRAPEProblem(
    t0 = total_time,
    ntimes = ntimes,
    ndim = ndim,
    H0 = H0,
    nb_additional_param = 0,
    error_sources = error_sources
)

# Define target unitary and subspace
target_unitary(x_add) = exp(-im * π/4 * (kron(σx, σx)))
projector = Matrix{Float64}(I, ndim, ndim)

fidelity_problem = FidelityRobustGRAPEProblem(
    unitary_problem = unitary_problem,
    projector = projector,
    target_unitary = target_unitary
)

# Create initial pulse amplitudes
x_initial = zeros(ntimes)  # Initial control parameters

# Define optimization parameters
fidelity_parameters = FidelityRobustGRAPEParameters(
    x_initial = x_initial,
    regularization_functions = [regularization_cost_deriv2],  # Smooth pulses
    regularization_coeff1 = [1e-5],  # Coefficient for first regularization
    regularization_coeff2 = [1e-5],  # Coefficient for second regularization
    error_source_coeff = [1e-3],     # Weight for error robustness
    iterations = 500,                # Maximum iterations
    solver_algorithm = LBFGS()       # Optimization algorithm
)

# Run optimization
result = optimize_fidelity_and_error_sources(fidelity_problem, fidelity_parameters)

# Extract optimized parameters
optimal_x = Optim.minimizer(result)

# Evaluate final fidelity
final_fidelity, _, _, _ = calculate_fidelity_and_derivatives(fidelity_problem, optimal_x)
println("Final fidelity: ", final_fidelity)
```

## Using the RydbergTools Module

RobustGRAPE.jl includes a specialized submodule for working with Rydberg atom systems, which provides pre-defined Hamiltonians and quantum gates commonly used in Rydberg atom quantum computing:

```julia
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra

# Access pre-defined Rydberg Hamiltonians
H = rydberg_hamiltonian_symmetric_blockaded(π/2, 0, 0)

# Create a CZ gate with a phase parameter
θ = π/4
cz_gate = cz_with_1q_phase_symmetric(θ)

# Define a Hamiltonian function for GRAPE optimization
function H0(t, ϕ, x_add)
    # ϕ contains the control parameters (laser phase)
    # x_add contains additional parameters
    return rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
end

# Now use this Hamiltonian in your GRAPE optimization
# ...
```

The RydbergTools module includes functions for:
- Hamiltonians in different Rydberg atom configurations (symmetric, blockaded, full)
- Standard quantum gates like the CZ gate with phase parameters
- Utility functions for phase handling (like unwrap_phase)

See the [Rydberg Tools API](../api/rydberg.md) for complete documentation of these functions.

## Workflow

The typical workflow with RobustGRAPE.jl involves:

1. **Define the quantum system**
   - Specify Hamiltonian and dimension
   - Define error sources
   - For Rydberg systems, use the RydbergTools module

2. **Set up the optimization problem**
   - Create problem definition
   - Specify target unitary operation
   - Configure regularization and robustness weights

3. **Run optimization**
   - Use `optimize_fidelity_and_error_sources` function
   - Adjust regularization parameters if needed

4. **Analyze results**
   - Calculate final fidelity
   - Analyze error robustness using `calculate_fidelity_response`
   - For frequency domain analysis, use `calculate_fidelity_response_fft`
   - Visualize optimized control pulses

## Example: Rydberg CZ Gate

Here's a simplified example of optimizing a CZ gate in a Rydberg atom system:

```julia
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Optim
using Random

# Define system parameters
ntimes = 200
total_time = 5.0

# Define Hamiltonian using RydbergTools
H0(t, ϕ, x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1], 0, 0)
cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])

# Create optimization problem
problem = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0 = total_time,
        ntimes = ntimes,
        ndim = 5,  # Dimension of the Rydberg system
        H0 = H0,
        nb_additional_param = 1,
        error_sources = []
    ),
    Diagonal([1, 2, 1, 0, 0]),  # Target state projection
    cz  # Target CZ operation
)

# Run optimization (simplified)
# ... 

# For complete examples, see the Examples section
```

## Next Steps

Explore the API documentation for more details on each component of RobustGRAPE.jl:

- [Types](../api/types.md): Core data structures
- [Unitary Calculations](../api/unitary.md): Computing evolution operators
- [Fidelity Calculations](../api/fidelity.md): Fidelity metrics and optimization
- [Regularization](../api/regularization.md): Pulse smoothing techniques
- [Rydberg Tools](../api/rydberg.md): Rydberg atom specific functionality

For complete working examples, check out the [Examples](../examples.md) page.