"""
    fidelity(ρ::AbstractMatrix, σ::AbstractMatrix)
    fidelity(ψ::AbstractVector, φ::AbstractVector)

Calculate the fidelity between two quantum states.

The states can be pure or mixed, and they are represented by vectors `ψ` and `φ` or density matrices `ρ` and `σ`, respectively.
"""
function fidelity(ρ::AbstractMatrix, σ::AbstractMatrix)
    abs2(tr(sqrt(ρ * σ)))
end

fidelity(ψ::AbstractVector, φ::AbstractVector) = abs2(ψ ⋅ φ)

function fidelity(ρ::AbstractArray, φ::AbstractVector)
    real(dot(φ, ρ, φ))
end

fidelity(ρ, φ) = fidelity(φ, ρ)

function get_w(λs, j)
    (sum(view(λs, 1:j)) - 1) / j
end

function project_onto_simplex!(λs)
    w = zero(eltype(λs))
    for i ∈ eachindex(λs)
        new_w = get_w(λs, i)
        if real(λs[i]) - new_w < 0
            break
        else
            w = new_w
        end
    end

    for (n, λ) in enumerate(λs)
        λs[n] = max(λ - w, zero(w))
    end

end

"""
    project2density(ρ)

Project a Hermitian matrix `ρ` to a density matrix by setting the negative eigenvalues to zero and normalizing the trace to 1.
"""
function project2density!(ρ)
    vals, vecs = eigen!(Hermitian(ρ), sortby=x -> -x)
    project_onto_simplex!(vals)
    broadcast!(√, vals, vals)

    rmul!(vecs, Diagonal(vals))
    mul!(ρ, vecs, vecs')
end

function project2density(ρ)
    σ = copy(ρ)
    project2density!(σ)
    σ
end

"""
    project2pure(ρ)

Project a Hermitian matrix `ρ` to a pure state by returning the eigenvector corresponding to the largest eigenvalue.
"""
function project2pure(ρ)
    F = eigen(Hermitian(ρ))
    F.vectors[:, end] # the last eigenvector is the one with the largest eigenvalue
end

"""
    polarization_state(s::Symbol, ::Type{T}=ComplexF32) where {S, T}

Return the polarization state corresponding to the symbol `S` as a vector of type `T`.
`S` can be one of the following symbols:
- `:H` for horizontal polarization
- `:V` for vertical polarization
- `:D` for diagonal polarization
- `:A` for antidiagonal polarization
- `:R` for right-handed circular polarization
- `:L` for left-handed circular polarization
"""
function polarization_state(s::Symbol, ::Type{T}=ComplexF32) where {T}
    if s == :H
        T[1, 0]
    elseif s == :V
        T[0, 1]
    elseif s == :D
        T[1/√2, 1/√2]
    elseif s == :A
        T[1/√2, -1/√2]
    elseif s == :R
        T[1/√2, -im/√2]
    elseif s == :L
        T[1/√2, im/√2]
    else
        throw(ArgumentError("Invalid polarization state symbol: $s. Expected one of :H, :V, :D, :A, :R, :L."))
    end
end

function get_noisy_probabilities(μ, θ, noise_level)
    probs = get_probabilities(μ, θ)
    for n ∈ eachindex(probs)
        shift = (2 * rand() - 1) * probs[n]
        probs[n] += shift * noise_level
    end
    normalize!(probs, 1)
    probs
end

"""
    isposdef!(ρ, xs, set)

Calculate the linear combination of the elements of `set` with the coefficients `xs` and check if the result is a positive definite matrix.
"""
function LinearAlgebra.isposdef!(ρ, θ)
    density_matrix_reconstruction!(ρ, θ)
    isposdef!(ρ)
end