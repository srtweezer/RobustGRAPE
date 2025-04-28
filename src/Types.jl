using Parameters
using Optim

"""
    ErrorSource

Represents a source of error in the Hamiltonian.

# Fields
- `Herror::Function`: A function that returns the error Hamiltonian. Signature: Herror(time, x, x_add, err)
"""
@with_kw struct ErrorSource
    Herror::Function # Hamiltonian(time,x,x_add,err)
end

"""
    UnitaryRobustGRAPEProblem

Represents a robust GRAPE optimization problem.

# Fields
- `t0::Real`: Total evolution time
- `ntimes::Int`: Number of time steps
- `ndim::Int`: Dimension of the Hilbert space
- `H0::Function`: Main Hamiltonian function. Signature: H0(time, x, x_add)
- `nb_additional_param::Int`: Number of additional parameters
- `error_sources::Vector{ErrorSource}`: List of error sources
- `ϵ::Real`: Small parameter for numerical differentiation
- `ϵ2::Real`: Second small parameter for numerical differentiation
"""
@with_kw struct UnitaryRobustGRAPEProblem
    t0::Real
    ntimes::Int
    ndim::Int
    H0::Function # Hamiltonian(time,x,x_add)
    nb_additional_param::Int
    error_sources::Vector{ErrorSource}
    ϵ::Real = 1e-8
    ϵ2::Real = 1e-4
end

"""
    FidelityRobustGRAPEProblem

Represents a robust GRAPE problem with fidelity calculation.

# Fields
- `unitary_problem::UnitaryRobustGRAPEProblem`: The underlying optimization problem
- `projector::Matrix{Real}`: Projector for subspace fidelity
- `target_unitary::Function`: Function that returns the target unitary. Signature: target_unitary(x_add)
"""
@with_kw struct FidelityRobustGRAPEProblem
    unitary_problem::UnitaryRobustGRAPEProblem
    projector::Matrix{Real}
    target_unitary::Function
end

"""
    FidelityRobustGRAPEParameters

Configuration parameters for robust quantum control optimization.

# Fields
- `x_initial::Vector{<:Real}`: Initial control pulse amplitudes and additional parameters
- `regularization_functions::Vector{Function}`: Functions to regularize the control pulses (one per parameter type)
- `regularization_coeff1::Vector{<:Real}`: First-order regularization coefficients
- `regularization_coeff2::Vector{<:Real}`: Second-order regularization coefficients
- `error_source_coeff::Vector{<:Real}`: Coefficients for each error source (must match error_sources in the problem)
- `time_limit::Real`: Maximum time for optimization in seconds (NaN means no limit)
- `iterations::Int`: Maximum number of optimization iterations
- `solver_algorithm::Optim.FirstOrderOptimizer`: Optimization algorithm (e.g., LBFGS(), GradientDescent())
- `additional_parameters::Dict{Symbol,Any}`: Additional parameters to pass to Optim.optimize
"""
@with_kw struct FidelityRobustGRAPEParameters
    x_initial::Vector{<:Real}
    regularization_functions::Vector{Function} # (nparam)
    regularization_coeff1::Vector{<:Real}
    regularization_coeff2::Vector{<:Real}
    error_source_coeff::Vector{<:Real}
    time_limit::Real = NaN
    iterations::Int = 1_000
    solver_algorithm::T where T<:Optim.FirstOrderOptimizer = Optim.LBFGS()
    additional_parameters::Dict{Symbol,Any} = Dict{Symbol,Any}()
end