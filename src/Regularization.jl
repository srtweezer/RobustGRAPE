"""
    regularization_cost(x::Vector{<:Real})

Compute regularization costs and their gradients to promote smoothness in control pulses.

This function calculates two regularization terms:
- First-order regularization (`reg1`): Penalizes large changes between consecutive elements (first derivative)
- Second-order regularization (`reg2`): Penalizes large changes in the rate of change (second derivative)

# Arguments
- `x::Vector{<:Real}`: The control parameter vector to regularize

# Returns
A tuple with four elements:
- `reg1::Real`: First-order regularization cost (sum of squared differences)
- `jac1::Vector{<:Real}`: Gradient of the first-order regularization with respect to `x`
- `reg2::Real`: Second-order regularization cost (sum of squared second differences)
- `jac2::Vector{<:Real}`: Gradient of the second-order regularization with respect to `x`

# Examples
```julia
x = [0.0, 0.1, 0.3, 0.2, 0.1]
reg1, jac1, reg2, jac2 = regularization_cost(x)
```
"""
function regularization_cost(x::Vector{<:Real})
    n = size(x,1)
    diff_x = diff(x)

    diff_diff_x = diff(diff_x)
    reg1 = sum(diff_x .^ 2)
    reg2 = sum(diff_diff_x .^ 2)

    jac1 = zeros(n)
    jac2 = zeros(n)
    jac1[2:n-1] .= -2 .* diff_diff_x
    jac1[1] += -2*diff_x[1]
    jac1[n] += 2*diff_x[n-1]
    jac2[1] = 2*(x[3]-2*x[2]+x[1])
    jac2[2] = 2*(x[4]-4*x[3]+5*x[2]-2*x[1])
    for i=3:n-2
        jac2[i] = 2*(x[i+2] - 4*x[i+1] + 6*x[i] - 4*x[i-1]+x[i-2])
    end
    jac2[n-1] = 2*(x[n-3]-4*x[n-2]+5*x[n-1]-2*x[n])
    jac2[n] = 2*(x[n-2]-2*x[n-1]+x[n])
    return reg1, jac1, reg2, jac2
end

"""
    regularization_cost(x::Vector{<:Real}, f::Function, df::Function)

Compute regularization costs and their gradients for transformed control parameters.

This function applies a transformation function `f` to the control parameters before
calculating regularization, and then applies the chain rule using `df` to compute
the gradient with respect to the original parameters.

# Arguments
- `x::Vector{<:Real}`: The original control parameter vector
- `f::Function`: Transformation function to apply to each element of `x`
- `df::Function`: Derivative of the transformation function

# Returns
A tuple with four elements:
- `reg1::Real`: First-order regularization cost on the transformed parameters
- `jac1::Vector{<:Real}`: Gradient of the first-order regularization with respect to original parameters
- `reg2::Real`: Second-order regularization cost on the transformed parameters
- `jac2::Vector{<:Real}`: Gradient of the second-order regularization with respect to original parameters

# Examples
```julia
x = [0.0, 0.1, 0.3, 0.2, 0.1]
f(x) = sin(x)
df(x) = cos(x)
reg1, jac1, reg2, jac2 = regularization_cost(x, f, df)
```
"""
function regularization_cost(x::Vector{<:Real}, f::Function, df::Function)
    f_x = f.(x)
    reg1, jac1, reg2, jac2 = regularization_cost(f_x)
    df_x = df.(x)
    return reg1, df_x .* jac1, reg2, df_x .* jac2
end

"""
    regularization_cost_phase(ϕs::Vector{<:Real})

Compute regularization costs and their gradients for phase-based control parameters.

This function calculates regularization terms for both sine and cosine of the phase values,
which promotes smoothness in the complex phasor representation of the phase controls.
This is particularly useful when optimizing phase-based control sequences where the physical
meaning depends on the periodic nature of phases.

# Arguments
- `ϕs::Vector{<:Real}`: The vector of phase values (in radians) to regularize

# Returns
A tuple with four elements:
- `reg1::Real`: First-order regularization cost (sum of cos and sin components)
- `jac1::Vector{<:Real}`: Gradient of the first-order regularization with respect to `ϕs`
- `reg2::Real`: Second-order regularization cost (sum of cos and sin components)
- `jac2::Vector{<:Real}`: Gradient of the second-order regularization with respect to `ϕs`

# Examples
```julia
ϕs = [0.0, 0.1, 0.3, 0.2, 0.1]
reg1, jac1, reg2, jac2 = regularization_cost_phase(ϕs)
```
"""
function regularization_cost_phase(ϕs::Vector{<:Real})
    reg1_cos, jac1_cos, reg2_cos, jac2_cos = regularization_cost(ϕs, x-> cos(x), x-> -sin(x))
    reg1_sin, jac1_sin, reg2_sin, jac2_sin = regularization_cost(ϕs, x-> sin(x), x-> cos(x))
    return reg1_cos+reg1_sin,jac1_cos+jac1_sin,reg2_cos+reg2_sin,jac2_cos+jac2_sin
end