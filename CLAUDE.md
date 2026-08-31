# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Package Overview

RobustGRAPE.jl is a Julia package for quantum optimal control with error robustness
analysis and optimization. It designs high-fidelity quantum gates that are robust against
noise and experimental imperfections, using GRAPE (GRadient Ascent Pulse Engineering).

- Optimizes pulse sequences with GRAPE
- Analyzes and optimizes robustness against error sources
- Computes fidelity response functions to characterize noise sensitivity
- Provides fidelity and robustness Hessians, and their principal control directions
- Ships Rydberg-specific Hamiltonians and CZ targets

## Interface

Hamiltonians and error Hamiltonians take a **discrete `time_step::Int`** between 1 and
`ntimes`, not a continuous time:

```julia
H0(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real})
Herror(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real}, err::Real)
```

They are passed as FunctionWrappers.jl wrappers, which keeps the inner propagation loop
type-stable: `HamiltonianFunctionWrapper`, `ErrorHamiltonianFunctionWrapper`,
`UnitaryFunctionWrapper`, `RegularizationFunctionWrapper`. Use the wrapped regularizers
(`regularization_cost_wrapped`, `regularization_cost_phase_wrapped`) rather than the bare
ones.

## Core Types

- `ErrorSource` — a source of error in the Hamiltonian
- `UnitaryRobustGRAPEProblem` — the control problem
- `FidelityRobustGRAPEProblem` — adds a target unitary and a projector
- `FidelityRobustGRAPEParameters` — optimization configuration

## Core Functions

- `calculate_unitary_and_derivatives` — evolution operator and its derivatives
- `calculate_fidelity_and_derivatives` — fidelity and error sensitivities
- `optimize_fidelity_and_error_sources` — high-level optimization
- `calculate_fidelity_response_fft` — fidelity response function to a noise PSD
- `calculate_fidelity_hessian`, `calculate_robustness_hessian`, `principal_waveforms` —
  low-rank structure of the landscape, for closed-loop calibration

## Examples

- Time-optimal CZ gate: `examples/time_optimal_cz.jl`
- Amplitude-robust CZ gate: `examples/ar_cz.jl`

## Build and test

```bash
julia --project=. -e "using Pkg; Pkg.test()"
julia --project=docs docs/make.jl        # docs; deployed from `stable` only
```

## Code style

- 4-space indentation
- Imports: Julia standard libraries, then external packages, then internal modules
- Type annotations on function arguments and struct fields where they help
- Docstrings with a Parameters section for public functions and types
- snake_case for functions and variables, CamelCase for types
- Export explicit symbols
- Descriptive error messages
