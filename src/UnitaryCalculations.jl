"""
    calculate_unitary_and_derivatives(problem::UnitaryRobustGRAPEProblem, x::Vector{<:Real})

Calculate the unitary evolution operator and its derivatives with respect to control parameters
and error sources.

# Arguments
- `problem::UnitaryRobustGRAPEProblem`: The robust GRAPE problem definition
- `x::Vector{<:Real}`: The optimization vector containing control parameters and additional parameters

# Returns
A tuple with:
- `U`: Final unitary evolution
- `U_dx`: Derivatives with respect to control parameters
- `U_dx_add`: Derivatives with respect to additional parameters
- `U_derr`: Derivatives with respect to error sources
- `U_derr_dx`: Mixed derivatives (error and control parameters)
- `U_derr_dx_add`: Mixed derivatives (error and additional parameters)
"""
function calculate_unitary_and_derivatives(problem::UnitaryRobustGRAPEProblem, x::Vector{<:Real})
    x_main = x[1:end-problem.nb_additional_param]
    @assert mod(size(x_main,1), problem.ntimes) == 0 "Control parameter size must be a multiple of time steps"
    ntimes = problem.ntimes
    nparam = size(x_main,1)÷ntimes
    x_main = reshape(x_main, (nparam,ntimes))
    x_add = x[end-problem.nb_additional_param+1:end]
    x_add_copy = copy(x_add)
    cum_evo = I
    old_cum_evo = copy(cum_evo)
    dt = problem.t0/problem.ntimes
    ndim = problem.ndim
    nerr = size(problem.error_sources,1)

    infimU_dx = zeros(Complex{typeof(x[1])},ndim,ndim,nparam,ntimes)
    infimU_dx_add = zeros(Complex{typeof(x[1])},ndim,ndim,problem.nb_additional_param,ntimes)
    infimU_derr = zeros(Complex{typeof(x[1])},ndim,ndim,nerr,ntimes)
    infimU_derr_dx = zeros(Complex{typeof(x[1])},ndim,ndim,nparam,nerr,ntimes)
    infimU_derr_dx_add = zeros(Complex{typeof(x[1])},ndim,ndim,problem.nb_additional_param,nerr,ntimes)

    infim_evo_derr_array = zeros(Complex{typeof(x[1])},ndim,ndim,nerr)
    infim_evo_dx_array = zeros(Complex{typeof(x[1])},ndim,ndim,nparam)
    infim_evo_dx_add_array = zeros(Complex{typeof(x[1])},ndim,ndim,problem.nb_additional_param)

    for nt=1:ntimes
        infim_evo = exp(-im*dt*problem.H0(nt, x_main[:,nt], x_add))
        cum_evo = infim_evo*cum_evo
        cum_evo_inv = inv(cum_evo)
        x_main_copy = copy(x_main[:,nt])
        for np=1:nparam
            x_main_copy[np] += problem.ϵ
            infim_evo_dx = exp(-im*dt*problem.H0(nt, x_main_copy, x_add))
            infimU_dx[:,:,np,nt] = cum_evo_inv*((1/problem.ϵ) *(infim_evo_dx-infim_evo))*old_cum_evo
            x_main_copy[np] = x_main[np,nt] + problem.ϵ2
            infim_evo_dx_array[:,:,np] = exp(-im*dt*problem.H0(nt, x_main_copy, x_add))
            x_main_copy[np] = x_main[np,nt]
        end
        for npa=1:problem.nb_additional_param
            x_add_copy[npa] += problem.ϵ
            infim_evo_dx_add = exp(-im*dt*problem.H0(nt, x_main[:,nt], x_add_copy))
            infimU_dx_add[:,:,npa,nt] = cum_evo_inv*((1/problem.ϵ) *(infim_evo_dx_add-infim_evo))*old_cum_evo
            x_add_copy[npa] = x_add[npa] + problem.ϵ2
            infim_evo_dx_add_array[:,:,npa] = exp(-im*dt*problem.H0(nt, x_main[:,nt], x_add_copy))
            x_add_copy[npa] = x_add[npa]
        end

        for ne=1:nerr
            infim_evo_derr = exp(-im*dt*(problem.error_sources[ne].Herror(nt,x_main[:,nt],x_add,problem.ϵ)
                + problem.H0(nt,x_main[:,nt],x_add))
            )
            infimU_derr[:,:,ne,nt] = cum_evo_inv*((1/problem.ϵ) * (infim_evo_derr-infim_evo))*old_cum_evo
            infim_evo_derr_array[:,:,ne] = exp(-im*dt*(problem.error_sources[ne].Herror(nt,x_main[:,nt],x_add,problem.ϵ2)
                + problem.H0(nt,x_main[:,nt],x_add))
            )

            for np=1:nparam
                x_main_copy[np] += problem.ϵ2
                infim_evo_derr_dx = exp(-im*dt*(problem.error_sources[ne].Herror(nt,x_main_copy,x_add,problem.ϵ2) +
                    problem.H0(nt,x_main_copy,x_add)
                ))
                infimU_derr_dx[:,:,np,ne,nt] = cum_evo_inv*((1/problem.ϵ2^2) * (
                    infim_evo_derr_dx + infim_evo
                    - infim_evo_derr_array[:,:,ne] - infim_evo_dx_array[:,:,np]
                ))*old_cum_evo
                x_main_copy[np] = x_main[np,nt]
            end

            for npa=1:problem.nb_additional_param
                x_add_copy[npa] += problem.ϵ2
                infim_evo_derr_dx_add = exp(-im*dt*(problem.error_sources[ne].Herror(nt,x_main[:,nt],x_add_copy,problem.ϵ2) +
                    problem.H0(nt,x_main[:,nt],x_add_copy)
                ))
                infimU_derr_dx_add[:,:,npa,ne,nt] = cum_evo_inv*((1/problem.ϵ2^2) * (
                    infim_evo_derr_dx_add + infim_evo
                    - infim_evo_derr_array[:,:,ne] - infim_evo_dx_add_array[:,:,npa]
                ))*old_cum_evo
                x_add_copy[npa] = x_add[npa]
            end
        end
        old_cum_evo = copy(cum_evo)
    end

    infimU_derr = permutedims(infimU_derr, (1,2,4,3))
    infimU_derr_dx = permutedims(infimU_derr_dx, (1,2,3,5,4))
    infimU_derr_dx_add = permutedims(infimU_derr_dx_add, (1,2,3,5,4))

    U_dx = zeros(Complex{typeof(x[1])},ndim,ndim,nparam,ntimes)
    U_dx_add = zeros(Complex{typeof(x[1])},ndim,ndim,problem.nb_additional_param)
    U_derr = zeros(Complex{typeof(x[1])},ndim,ndim,nerr)
    U_derr_dx = zeros(Complex{typeof(x[1])},ndim,ndim,nparam,ntimes,nerr)
    U_derr_dx_add = zeros(Complex{typeof(x[1])},ndim,ndim,problem.nb_additional_param,nerr)

    infimU_derr_cumsum = cumsum(infimU_derr,dims=3)
    infimU_derr_revcumsum = reverse(cumsum(reverse(infimU_derr,dims=3),dims=3),dims=3)
    for nt=1:ntimes
        for np=1:nparam
            U_dx[:,:,np,nt] = cum_evo*infimU_dx[:,:,np,nt]
        end
    end
    for npa=1:problem.nb_additional_param
        U_dx_add[:,:,npa] = cum_evo*sum(infimU_dx_add[:,:,npa,:],dims=3)[:,:,1]
    end
    for ne=1:nerr
        U_derr[:,:,ne] = cum_evo*sum(infimU_derr[:,:,:,ne],dims=3)[:,:,1]
        for nt=2:ntimes
            for np=1:nparam
                U_derr_dx[:,:,np,nt,ne] += infimU_dx[:,:,np,nt]*infimU_derr_cumsum[:,:,nt-1,ne]
            end
        end
        for nt=1:(ntimes-1)
            for np=1:nparam
                U_derr_dx[:,:,np,nt,ne] += infimU_derr_revcumsum[:,:,nt+1,ne]*infimU_dx[:,:,np,nt]
            end
        end
        for nt=1:ntimes
            for np=1:nparam
                U_derr_dx[:,:,np,nt,ne] += infimU_derr_dx[:,:,np,nt,ne]
                U_derr_dx[:,:,np,nt,ne] = cum_evo*U_derr_dx[:,:,np,nt,ne]
            end
        end
        for npa=1:problem.nb_additional_param
            for nt=2:ntimes
                U_derr_dx_add[:,:,npa,ne] += infimU_dx_add[:,:,npa,nt]*infimU_derr_cumsum[:,:,nt-1,ne]
            end
            for nt=1:(ntimes-1)
                U_derr_dx_add[:,:,npa,ne] += infimU_derr_revcumsum[:,:,nt+1,ne]*infimU_dx_add[:,:,npa,nt]
            end
            for nt=1:ntimes
                U_derr_dx_add[:,:,npa,ne] += infimU_derr_dx_add[:,:,npa,nt,ne]
            end
            U_derr_dx_add[:,:,npa,ne] = cum_evo*U_derr_dx_add[:,:,npa,ne]
        end
    end

    return (cum_evo, U_dx, U_dx_add, U_derr, U_derr_dx, U_derr_dx_add)
end

"""
    calculate_interaction_error_operators(problem::UnitaryRobustGRAPEProblem, x::Vector{<:Real})

Calculate the interaction picture representation of error operators at each time step.

This function transforms error operators from the Schrödinger picture to the interaction picture,
which is essential for analyzing how errors affect the quantum dynamics throughout the evolution.
The interaction picture provides a way to separate the influence of the control Hamiltonian from
the error terms.

# Parameters
- `problem::UnitaryRobustGRAPEProblem`: The robust GRAPE problem definition
- `x::Vector{<:Real}`: The optimization vector containing control parameters and additional parameters

# Returns
- A tensor of dimensions (ndim, ndim, ntimes, nerr) containing the interaction picture
  representation of each error operator at each time step.

# Notes
- The interaction picture transformation uses the cumulative evolution operator
- Error operators are scaled by the small parameter ϵ used for numerical differentiation
- The returned tensor has dimensions permuted for convenient access to time-dependent error operators
"""
function calculate_interaction_error_operators(problem::UnitaryRobustGRAPEProblem, x::Vector{<:Real})
    x_main = x[1:end-problem.nb_additional_param]
    @assert mod(size(x_main,1), problem.ntimes) == 0 "Control parameter size must be a multiple of time steps"
    ntimes = problem.ntimes
    nparam = size(x_main,1)÷ntimes
    x_main = reshape(x_main, (nparam,ntimes))
    x_add = x[end-problem.nb_additional_param+1:end]
    cum_evo = I
    dt = problem.t0/problem.ntimes
    ndim = problem.ndim
    nerr = size(problem.error_sources,1)

    error_operators_int = zeros(Complex{typeof(x[1])},ndim,ndim,nerr,ntimes)
    for nt=1:ntimes
        cum_evo_inv = inv(cum_evo)
        for ne=1:nerr
            Oerr = (1/problem.ϵ) * problem.error_sources[ne].Herror(nt,x_main[:,nt],x_add,problem.ϵ)
            error_operators_int[:,:,ne,nt] = cum_evo_inv*Oerr*cum_evo
        end
        infim_evo = exp(-im*dt*problem.H0(nt, x_main[:,nt], x_add))
        cum_evo = infim_evo*cum_evo
    end

    return permutedims(error_operators_int, (1,2,4,3))
end
