# Types

```@meta
CurrentModule = RobustGRAPE
```

## Core Types

```@docs
ErrorSource
UnitaryRobustGRAPEProblem
FidelityRobustGRAPEProblem
FidelityRobustGRAPEParameters
```
## Function Wrappers

Hamiltonians, error Hamiltonians, target unitaries and regularizers are passed as
[FunctionWrappers.jl](https://github.com/yuyichao/FunctionWrappers.jl) wrappers rather
than as bare functions, which keeps the inner propagation loop type-stable. Note that
the Hamiltonian and error signatures take a **`time_step::Int`** between 1 and `ntimes`,
not a continuous time:

```julia
H0(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real})
Herror(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real}, err::Real)
```

Use the wrapped regularizers (`regularization_cost_wrapped`,
`regularization_cost_phase_wrapped`) in the same way.

```@docs
HamiltonianFunctionWrapper
ErrorHamiltonianFunctionWrapper
UnitaryFunctionWrapper
RegularizationFunctionWrapper
```
