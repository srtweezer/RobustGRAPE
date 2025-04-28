# RobustGRAPE.jl

A Julia package for robust quantum optimal control using GRadient Ascent Pulse Engineering (GRAPE).

## Overview

RobustGRAPE.jl provides tools for designing high-fidelity quantum gates that are robust against noise and experimental errors. It extends the standard GRAPE algorithm to optimize quantum control pulses while minimizing their sensitivity to various error sources.

The package is particularly suited for optimizing quantum gates in systems such as:
- Rydberg atom arrays
- Trapped ions
- Superconducting qubits
- NMR systems

## Features

- Optimization of quantum control pulses for high-fidelity gates
- Robustness against various error sources and noise
- Frequency domain analysis of pulse sensitivity
- Support for arbitrary Hamiltonians and quantum systems
- Efficient gradient calculations for optimization
- Various regularization techniques for pulse smoothing

## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add("RobustGRAPE")
```

Or to install the development version from this repository:

```julia
using Pkg
Pkg.add(url="https://github.com/USERNAME/RobustGRAPE.jl")
```

## Quick Start

```julia
using RobustGRAPE
using RobustGRAPE.RydbergTools
using LinearAlgebra
using Optim

# Define a quantum system and target operation
H0(t,ϕ,x_add) = rydberg_hamiltonian_symmetric_blockaded(ϕ[1],0,0)
cz(x_add) = cz_with_1q_phase_symmetric(x_add[1])

# Create a GRAPE optimization problem
problem = FidelityRobustGRAPEProblem(
    UnitaryRobustGRAPEProblem(
        t0=7.0,             # Total evolution time
        ntimes=500,         # Number of timesteps
        ndim=5,             # Dimension of the system
        H0=H0,              # Hamiltonian function
        nb_additional_param=1,
        error_sources=[]    # No error sources for basic optimization
    ),
    Diagonal([1,2,1,0,0]),  # Target state projection
    cz                      # Target operation
)

# Configure optimization parameters
params = FidelityRobustGRAPEParameters(
    x_initial = rand(501),  # 500 timesteps + 1 phase parameter
    regularization_functions = [regularization_cost_phase],
    regularization_coeff1=[1e-7],
    regularization_coeff2=[1e-7],
    time_limit=40
)

# Run optimization
result = optimize_fidelity_and_error_sources(problem, params)
```

## Documentation

For more information, see the [online documentation](https://srtweezer.github.io/RobustGRAPE.jl/stable/).

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This package was developed with assistance from [Claude](https://claude.ai), particularly for documentation and code organization.