using Parameters
using Optim

"""
    ErrorSource

Represents a source of error in the Hamiltonian.

# Fields
- `Herror::Function`: A function that returns the error Hamiltonian matrix. Signature: `Herror(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real}, err::Real)`. `x` has for size the number of main parameters, `x_add` has for size the number of additional parameters. `time_step` is an integer between 1 and ntimes.
"""
@with_kw struct ErrorSource
    Herror::Function
end

"""
    UnitaryRobustGRAPEProblem

Represents a robust GRAPE unitary calculation problem.

# Fields
- `t0::Real`: Total evolution time
- `ntimes::Int`: Number of time steps
- `ndim::Int`: Dimension of the Hilbert space
- `H0::Function`: Main Hamiltonian function (must return a matrix). Signature: `H0(time_step::Int, x::Vector{<:Real}, x_add::Vector{<:Real})`. `x` has for size the number of main parameters, `x_add` has for size the number of additional parameters. `time_step` is an integer between 1 and ntimes.
- `nb_additional_param::Int`: Number of additional parameters
- `error_sources::Vector{ErrorSource}`: List of error sources
- `ϵ::Real = 1e-8`: Small parameter for first-order finite difference
- `ϵ2::Real = 1e-4`: Small parameter for second-order finite difference
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

Represents a robust GRAPE fidelity calculation problem.

# Fields
- `unitary_problem::UnitaryRobustGRAPEProblem`: The underlying optimization problem
- `projector::Matrix{<:Real}`: Projector for subspace fidelity. It can be a "pseudo"-projector (see Quick start and example for more information)
- `target_unitary::Function`: Function that returns the target unitary. Signature: `target_unitary(x_add::Vector{<:Real})` where `x_add` has size the number of additional parameters.
"""
@with_kw struct FidelityRobustGRAPEProblem
    unitary_problem::UnitaryRobustGRAPEProblem
    projector::Matrix{<:Real}
    target_unitary::Function
end

"""
    FidelityRobustGRAPEParameters

Configuration parameters for parameterized quantum gate optimization according to the GRAPE algorithm with optional robustness.

# Fields
- `x_initial::Vector{<:Real}`: Initial control pulse amplitudes and additional parameters
- `regularization_functions::Vector{Function}`: Vector of functions to regularize the control pulses. The vector has size the number of main parameters. Each function has signature `(x::Vector{<:Real}) -> (r1::Real,j1::Vector{<:Real},r2::Real,j2::Vector{<:Real})` where `r1` and `r2` are the first and second-order costs, and `j1` and `j2` are the gradients of these costs w.r.t. the control parameters. `x`, `j1` and `j2` have size `ntimes`.
- `regularization_coeff1::Vector{<:Real}`: First-order regularization coefficients; has size the number of main parameters.
- `regularization_coeff2::Vector{<:Real}`: Second-order regularization coefficients; has size the number of main parameters.
- `error_source_coeff::Vector{<:Real}`: Coefficients for each error source. Has size the number of error sources.
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