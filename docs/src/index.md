# RobustGRAPE.jl

*A Julia package for robust quantum optimal control using the GRAPE algorithm*

## Overview

RobustGRAPE.jl is a package for optimizing quantum control pulses using the GRadient Ascent Pulse Engineering (GRAPE) algorithm, with additional features for robustness against noise and errors. The package is designed for optimizing control pulses for quantum gates in systems like Rydberg atoms, trapped ions, and superconducting qubits.

## Features

- Optimization of quantum control pulses for high-fidelity gates
- Robustness against various error sources and noise
- Frequency domain analysis of robustness
- Support for arbitrary Hamiltonians and quantum systems
- Efficient gradient calculations for optimization
- Various regularization techniques for pulse smoothing
- Specialized tools for Rydberg atom systems through the RydbergTools module

## RydbergTools Module

The package includes a specialized module for working with Rydberg atom systems:

```julia
using RobustGRAPE.RydbergTools

# Access pre-defined Rydberg Hamiltonians
H = rydberg_hamiltonian_symmetric_blockaded(π/2, 0, 0)

# Create a CZ gate with a phase parameter
cz_gate = cz_with_1q_phase_symmetric(π/4)
```

For more details, see the [Rydberg Tools API](api/rydberg.md).

## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add("RobustGRAPE")
```

Or to install the development version from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/USERNAME/RobustGRAPE.jl")
```

## Getting Started

See the [Getting Started](guide/getting-started.md) guide for a quick introduction to using RobustGRAPE.jl.

## Examples

Check out the [Examples](examples.md) page for complete working examples, including:

- Time-optimal CZ gate implementation using Rydberg atoms
- Robust pulse design against intensity and frequency errors
- Frequency domain analysis of optimized pulses

## Acknowledgments

This package's documentation and code organization was developed with the assistance of [Claude](https://claude.ai), an AI assistant by Anthropic.