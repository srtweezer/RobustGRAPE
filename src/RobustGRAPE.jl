module RobustGRAPE

using LinearAlgebra
using FunctionWrappers: FunctionWrapper

# Export main types
export ErrorSource, UnitaryRobustGRAPEProblem, FidelityRobustGRAPEProblem, FidelityRobustGRAPEParameters
export HamiltonianFunctionWrapper, ErrorHamiltonianFunctionWrapper, RegularizationFunctionWrapper, UnitaryFunctionWrapper

# Export main functions
export calculate_unitary_and_derivatives, calculate_fidelity_and_derivatives, optimize_fidelity_and_error_sources
export calculate_interaction_error_operators, calculate_fidelity_response, calculate_fidelity_response_fft, calculate_expectation_values

# Export regularization
export regularization_cost, regularization_cost_phase
export regularization_cost_wrapped, regularization_cost_phase_wrapped

include("Types.jl")
include("UnitaryCalculations.jl")
include("FidelityCalculations.jl")
include("Regularization.jl")

module RydbergTools
    include("RydbergTools.jl")
end


end # module RobustGRAPE
