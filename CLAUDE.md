# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview
RobustGRAPE.jl is a Julia package for quantum optimal control with error robustness analysis and optimization. It focuses on designing high-fidelity quantum gates that are robust against noise and experimental imperfections using the GRAPE (GRadient Ascent Pulse Engineering) algorithm.

Key features:
- Optimizes quantum gate pulse sequences using GRAPE algorithm
- Analyzes and optimizes robustness against various error sources
- Calculates fidelity response functions to characterize sensitivity to noise
- Provides tools specifically for Rydberg atom-based quantum gates

## Important Interface Changes
- **Hamiltonian and error functions now use `time_step::Int` instead of `time::Real`**
- Function signatures are:
  - `H0(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real})`
  - `Herror(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real}, err::Real)`

## Core Types
- `ErrorSource`: Represents a source of error in the Hamiltonian
- `UnitaryRobustGRAPEProblem`: Defines the quantum control problem parameters
- `FidelityRobustGRAPEProblem`: Adds target unitary and projector to calculate fidelity
- `FidelityRobustGRAPEParameters`: Configuration for optimization (regularization, error coefficients, etc.)

## Core Functions 
- `calculate_unitary_and_derivatives`: Computes the evolution operator and its derivatives
- `calculate_fidelity_and_derivatives`: Computes fidelity and sensitivity to errors
- `optimize_fidelity_and_error_sources`: High-level optimization interface
- `calculate_fidelity_response_fft`: Computes the fidelity response function to noise

## Examples
The package includes two key examples:
- Time-optimal CZ gate (`examples/time_optimal_cz.jl`)
- Amplitude-robust CZ gate (`examples/ar_cz.jl`)

## Build/Test Commands
- Run all tests: `julia --project=. -e "using Pkg; Pkg.test()"`
- Run single test: `julia --project=. -e "using Pkg; Pkg.test(\"RobustGRAPE\", test_args=[\"path/to/test.jl\"])"`

## Code Style Guidelines
- Use 4-space indentation
- Import order: Julia standard libraries first, then external packages, then internal modules
- Type annotations should be used for function arguments and struct fields when helpful
- Use docstrings with Parameters section for public functions and types
- Error handling: use descriptive error messages
- Follow Julia naming conventions: snake_case for functions and variables, CamelCase for types
- Export explicit symbols rather than `export *` patterns
- Parameters.jl `@with_kw` is used selectively, not required for all structs