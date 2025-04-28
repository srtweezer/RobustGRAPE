# Getting Started

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

## Workflow

The typical workflow with RobustGRAPE.jl involves:

1. **Define the quantum system**
   - Specify Hamiltonian and dimension
   - Define error sources

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
   - Visualize optimized control pulses

## Next Steps

Explore the API documentation for more details on each component of RobustGRAPE.jl:

- [Types](../api/types.md): Core data structures
- [Unitary Calculations](../api/unitary.md): Computing evolution operators
- [Fidelity Calculations](../api/fidelity.md): Fidelity metrics and optimization
- [Regularization](../api/regularization.md): Pulse smoothing techniques
- [Rydberg Tools](../api/rydberg.md): Rydberg atom specific functionality