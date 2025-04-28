# RobustGRAPE Tests

This directory contains tests for the RobustGRAPE package.

## Running Tests

To run all tests, use Julia's built-in test system:

```julia
julia> using Pkg
julia> Pkg.test("RobustGRAPE")
```

Or, from the RobustGRAPE directory:

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

## Test Descriptions

### Time-optimal CZ Gate Test

This test verifies the functionality of the robust GRAPE implementation using a known solution for a time-optimal CZ gate with Rydberg atoms.

#### Physical Background

The test implements a quantum control protocol for generating a controlled-Z (CZ) gate between two Rydberg atoms. The control is achieved by varying the phase of the driving field over time, following a specific pattern:

φ(t) = A·cos(ω₀·t - φ₀) + δ₀·t

where:
- A: Amplitude of the oscillation
- ω₀: Frequency of the oscillation
- φ₀: Phase offset
- δ₀: Linear detuning slope

This control shape, combined with a specific single-qubit phase correction (θ), is known to create a high-fidelity CZ gate in the computational subspace of the two-atom system.

#### Test Implementation

The test:
1. Sets up a control pulse with known optimal parameters
2. Defines the Hamiltonian using `rydberg_hamiltonian_symmetric_blockaded`
3. Creates a `FidelityRobustGRAPEProblem` instance
4. Calculates the gate fidelity using `calculate_fidelity_and_derivatives`
5. Verifies the fidelity exceeds 0.9999 (infidelity < 10⁻⁴)

This test ensures that the core fidelity calculation and unitary evolution functions are working correctly.

#### Parameters

The specific parameter values used are:
- t₀ = 2π·1.22 (total gate time)
- A = 0.7701624 (amplitude)
- ω₀ = 0.97525275 (frequency)
- φ₀ = -0.97449603 (phase offset)
- δ₀ = -0.04319765 (detuning slope)
- θ = 2.0802725844516097 (single-qubit phase correction)

These parameters were found and described in the following reference:
Evered, S.J., Bluvstein, D., Kalinowski, M. et al. High-fidelity parallel entangling gates on a neutral-atom quantum computer. Nature 622, 268–272 (2023).

### Fidelity Gradient Validation Test

This test validates the correctness of the analytical gradient calculation in the `calculate_fidelity_and_derivatives` function using finite difference approximation.

The test:
1. Sets up a Rydberg atom quantum control problem
2. Generates random control parameters
3. Compares the analytical gradient from the package with a numerical approximation

By testing both control parameters and the additional phase parameter, this test ensures that the gradient calculations are correctly implemented, which is essential for optimization algorithms used in quantum optimal control.

### Gradient-based Pulse Optimization Test

This test implements gradient-based optimization to find a high-fidelity pulse for a CZ gate using the LBFGS optimization algorithm from the Optim.jl package.

#### Implementation

The test demonstrates the use of RobustGRAPE for quantum optimal control:

1. Sets up a reduced-size problem (fewer time steps for efficiency)
2. Implements the cost function that combines:
   - Infidelity (1 - fidelity) as the primary objective
   - Phase regularization to encourage smooth pulses
3. Runs LBFGS optimization for up to 200 iterations
4. Verifies that the optimizer finds a high-fidelity solution (infidelity < 10⁻⁶)

This test validates RobustGRAPE's ability to find optimal control pulses from a random starting point, demonstrating its practical utility for quantum gate design. The implementation closely mirrors the approach in the time_optimal_cz.ipynb example notebook.