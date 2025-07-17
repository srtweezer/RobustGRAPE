using LinearAlgebra

"""
    rydberg_hamiltonian_symmetric_blockaded(ϕ::Real, ϵ::Real, δ::Real)

Constructs the Hamiltonian for a symmetric Rydberg-blockaded two-atom system.

# Basis
`|00⟩, |01⟩, |11⟩, |0r⟩, |W⟩` where `|W⟩ = (|1r⟩ + |r1⟩)/√2`

# Parameters
- `ϕ::Real`: Phase of the driving field
- `ϵ::Real`: Relative amplitude deviation parameter
- `δ::Real`: Detuning of the Rydberg state

# Mathematical form
```math
H = 
\\begin{pmatrix}
0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & \\frac{(1+\\epsilon)e^{-i\\phi}}{2} & 0 \\\\
0 & 0 & 0 & 0 & \\frac{(1+\\epsilon)e^{-i\\phi}}{\\sqrt{2}} \\\\
0 & \\frac{(1+\\epsilon)e^{i\\phi}}{2} & 0 & \\delta & 0 \\\\
0 & 0 & \\frac{(1+\\epsilon)e^{i\\phi}}{\\sqrt{2}} & 0 & \\delta
\\end{pmatrix}
```

# Returns
- Matrix representing the Hamiltonian in the symmetric basis described above
"""
function rydberg_hamiltonian_symmetric_blockaded(ϕ::Real,ϵ::Real,δ::Real)
    return [
        0 0 0 0 0
        0 0 0 exp(-im*ϕ)*(1+ϵ)/2 0
        0 0 0 0 exp(-im*ϕ)*(1+ϵ)/√2
        0 exp(im*ϕ)*(1+ϵ)/2 0 δ 0
        0 0 exp(im*ϕ)*(1+ϵ)/√2 0 δ
    ]
end

"""
    rydberg_hamiltonian_full_blockaded(ϕ::Real, ϵ::Real, δ::Real)

Constructs the Hamiltonian for a fully-described Rydberg-blockaded two-atom system.

# Basis
`|00⟩, |01⟩, |10⟩, |11⟩, |0r⟩, |r0⟩, |W'⟩` where `|W'⟩ = (|1r⟩ + |r1⟩)/√2`

# Parameters
- `ϕ::Real`: Phase of the driving field
- `ϵ::Real`: Relative amplitude deviation parameter
- `δ::Real`: Detuning of the Rydberg state

# Mathematical form
```math
H = 
\\begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & \\frac{(1+\\epsilon)e^{-i\\phi}}{2} & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & \\frac{(1+\\epsilon)e^{-i\\phi}}{2} & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & \\frac{(1+\\epsilon)e^{-i\\phi}}{\\sqrt{2}} \\\\
0 & \\frac{(1+\\epsilon)e^{i\\phi}}{2} & 0 & 0 & \\delta & 0 & 0 \\\\
0 & 0 & \\frac{(1+\\epsilon)e^{i\\phi}}{2} & 0 & 0 & \\delta & 0 \\\\
0 & 0 & 0 & \\frac{(1+\\epsilon)e^{i\\phi}}{\\sqrt{2}} & 0 & 0 & \\delta
\\end{pmatrix}
```

# Returns
- Matrix representing the Hamiltonian in the basis described above
"""
function rydberg_hamiltonian_full_blockaded(ϕ::Real,ϵ::Real,δ::Real)
    return [
        0 0 0 0 0 0 0
        0 0 0 0 exp(-im*ϕ)*(1+ϵ)/2 0 0
        0 0 0 0 0 exp(-im*ϕ)*(1+ϵ)/2 0
        0 0 0 0 0 0 exp(-im*ϕ)*(1+ϵ)/√2
        0 exp(im*ϕ)*(1+ϵ)/2 0 0 δ 0 0
        0 0 exp(im*ϕ)*(1+ϵ)/2 0 0 δ 0
        0 0 0 exp(im*ϕ)*(1+ϵ)/√2 0 0 δ
    ]
end

"""
    rydberg_hamiltonian_full(ϕ::Real, Ω1::Real, Ω2::Real, δ1::Real, δ2::Real, B::Real)

Constructs the full Hamiltonian for a two-atom Rydberg system without symmetry constraints.

# Basis
`|00⟩, |01⟩, |10⟩, |11⟩, |0r⟩, |r0⟩, |1r⟩, |r1⟩, |rr⟩`

# Parameters
- `ϕ::Real`: Phase of the driving field
- `Ω1::Real`: Rabi frequency for the first atom
- `Ω2::Real`: Rabi frequency for the second atom
- `δ1::Real`: Detuning for the first atom
- `δ2::Real`: Detuning for the second atom
- `B::Real`: Rydberg-Rydberg blockade shift

# Mathematical form
```math
H = 
\\begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & \\frac{\\Omega_1 e^{-i\\phi}}{2} & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & \\frac{\\Omega_2 e^{-i\\phi}}{2} & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & \\frac{\\Omega_1 e^{-i\\phi}}{2} & \\frac{\\Omega_2 e^{-i\\phi}}{2} & 0 \\\\
0 & \\frac{\\Omega_1 e^{i\\phi}}{2} & 0 & 0 & \\delta_1 & 0 & 0 & 0 & 0 \\\\
0 & 0 & \\frac{\\Omega_2 e^{i\\phi}}{2} & 0 & 0 & \\delta_2 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & \\frac{\\Omega_1 e^{i\\phi}}{2} & 0 & \\delta_1 & 0 & \\frac{\\Omega_2 e^{-i\\phi}}{2} \\\\
0 & 0 & 0 & 0 & 0 & \\frac{\\Omega_2 e^{i\\phi}}{2} & 0 & \\delta_2 & \\frac{\\Omega_1 e^{-i\\phi}}{2} \\\\
0 & 0 & 0 & 0 & 0 & 0 & \\frac{\\Omega_2 e^{i\\phi}}{2} & \\frac{\\Omega_1 e^{i\\phi}}{2} & \\delta_1 + \\delta_2 + B
\\end{pmatrix}
```

# Returns
- Matrix representing the Hamiltonian in the full basis described above
"""
function rydberg_hamiltonian_full(ϕ::Real,Ω1::Real,Ω2::Real,δ1::Real,δ2::Real,B::Real)
    return [
        0 0 0 0 0 0 0 0 0
        0 0 0 0 exp(-im*ϕ)*Ω1/2 0 0 0 0
        0 0 0 0 0 exp(-im*ϕ)*Ω2/2 0 0 0
        0 0 0 0 0 0 exp(-im*ϕ)*Ω1/2 exp(-im*ϕ)*Ω2/2 0
        0 exp(im*ϕ)*Ω1/2 0 0 δ1 0 0 0 0
        0 0 exp(im*ϕ)*Ω2/2 0 0 δ2 0 0 0
        0 0 0 exp(im*ϕ)*Ω1/2 0 0 δ1 0 exp(-im*ϕ)*Ω2/2
        0 0 0 exp(im*ϕ)*Ω2/2 0 0 0 δ2 exp(-im*ϕ)*Ω1/2
        0 0 0 0 0 0 exp(im*ϕ)*Ω2/2 exp(im*ϕ)*Ω1/2 δ1+δ2+B
    ]
end

"""    
    cz_with_1q_phase_symmetric(θ::Real)

Constructs the CZ gate with additional single-qubit phase in the symmetric subspace.

# Basis
`|00⟩, |01⟩, |11⟩, |0r⟩, |W⟩` (same as in `rydberg_hamiltonian_symmetric_blockaded`)

# Parameters
- `θ::Real`: Single-qubit phase parameter

# Mathematical form
```math
U = 
\\begin{pmatrix}
1 & 0 & 0 & 0 & 0 \\\\
0 & e^{i\\theta} & 0 & 0 & 0 \\\\
0 & 0 & e^{i(2\\theta+\\pi)} & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0
\\end{pmatrix}
```

The diagonal structure encodes a CZ gate with additional single-qubit phase rotations.

# Returns
- Diagonal matrix representing the CZ gate with phase rotations
"""
function cz_with_1q_phase_symmetric(θ::Real)
    return collect(ComplexF64,Diagonal([1,exp(im*θ),exp(im*(2*θ+π)),0,0]))
end

"""    
    cz_with_1q_phase_full(θ::Real; rydberg_dimension::Int = 5)

Constructs the CZ gate with additional single-qubit phase in the full computational basis.

# Basis
`|00⟩, |01⟩, |10⟩, |11⟩, |0r⟩, |r0⟩, |1r⟩, |r1⟩, |rr⟩` (same as in `rydberg_hamiltonian_full`)

# Parameters
- `θ::Real`: Single-qubit phase parameter
- `rydberg_dimension::Int=5`: Dimension of the Rydberg subspace (optional, default: 5)

# Mathematical form
```math
U = 
\\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & e^{i\\theta} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & e^{i\\theta} & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & e^{i(2\\theta+\\pi)} & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\\end{pmatrix}
```

The diagonal structure encodes a CZ gate with additional single-qubit phase rotations.

# Returns
- Diagonal matrix representing the CZ gate with phase rotations in the full basis
"""
function cz_with_1q_phase_full(θ::Real; rydberg_dimension::Int = 5)
    diag_vec = zeros(Complex{typeof(θ)},4+rydberg_dimension)
    diag_vec[1] = 1
    diag_vec[2:3] .= exp(im*θ)
    diag_vec[4] = exp(im*(2*θ+π))
    return collect(Complex{typeof(θ)},Diagonal(diag_vec))
end


"""
    unwrap_phase(ϕ)

Unwraps a sequence of phase values by removing jumps greater than π.

Adjusts phase values to maintain continuity across the 2π boundary, eliminating
artificial discontinuities in phase data while preserving the actual phase evolution.
Particularly useful for plotting phase values to avoid discontinuous jumps.

# Parameters
- `ϕ`: Array of phase values to unwrap

# Returns
- An array of unwrapped phase values with the same length as the input
"""
function unwrap_phase(ϕ)
    ϕ2 = copy(ϕ)
    ϕ2 = mod.(ϕ2, 2*π)
    for i=1:(size(ϕ2,1)-1)
        if ϕ2[i+1]-ϕ2[i] > π
            ϕ2[i+1:end] .-= 2*π
        elseif ϕ2[i+1]-ϕ2[i] < -π
            ϕ2[i+1:end] .+= 2*π
        end
    end
    return ϕ2
end

export rydberg_hamiltonian_symmetric_blockaded
export rydberg_hamiltonian_full_blockaded
export rydberg_hamiltonian_full
export cz_with_1q_phase_symmetric
export cz_with_1q_phase_full
export unwrap_phase