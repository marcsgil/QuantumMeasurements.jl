abstract type QuantumMeasurementRandom end

"""
    HaarUnitary(dim::Int)

A type representing a Haar-random unitary matrix of dimension `dim`.
"""
struct HaarUnitary <: QuantumMeasurementRandom
    dim::Int
end

"""
    Base.rand([rng = default_rng()], type, [n_samples])

Sample `n_samples` from `type`.

If `n_samples` is not provided, a single sample is returned.

Possible values for type are [`HaarUnitary`](@ref), [`HaarVector`](@ref), [`Simplex`](@ref), [`ProductMeasure`](@ref), and [`GinibreEnsamble`](@ref).
"""
function Base.rand(rng::AbstractRNG, type::HaarUnitary, n_samples::Integer)
    Zs = randn(rng, ComplexF32, type.dim, type.dim, n_samples)

    for Z ∈ eachslice(Zs, dims=3)
        Q, R = qr(Z)
        Λ = diag(R)
        @. Λ /= abs(Λ)
        Z .= Q * diagm(Λ)
    end

    return Zs
end

"""
    HaarVector(dim::Int)

A type representing a Haar-random unitary vector of dimension `dim`.
"""
struct HaarVector  <: QuantumMeasurementRandom
    dim::Int
end

function Base.rand(rng::AbstractRNG, type::HaarVector, n_samples::Integer)
    rand(rng, HaarUnitary(type.dim), n_samples)[1, :, :]
end

"""
    Simplex(dim::Int)

A type representing a random point on the simplex embeded in a space of dimension `dim`.
"""
struct Simplex  <: QuantumMeasurementRandom
    dim::Int
end

function Base.rand(rng::AbstractRNG, type::Simplex)
    dim = type.dim
    ξs = rand(rng, Float32, dim - 1)
    λs = Vector{Float32}(undef, dim)
    for (k, ξ) ∈ enumerate(ξs)
        λs[k] = (1 - ξ^(1 / (dim - k))) * (1 - sum(λs[1:k-1]))
    end
    λs[end] = 1 - sum(λs[1:end-1])
    λs
end

function Base.rand(rng::AbstractRNG, type::Simplex, nsamples::Integer)
    stack(rand(rng, type) for _ ∈ 1:nsamples)
end

function combine!(unitaries::Array{T1,3}, probabilities::Array{T2,2}) where {T1,T2}
    @assert size(probabilities, 2) == size(unitaries, 3) "The number of probabilities and unitaries must be the same."
    @assert size(probabilities, 1) == size(unitaries, 1) "The dimension of the probabilities and unitaries must be the same."

    for (U, p) ∈ zip(eachslice(unitaries, dims=3), eachslice(probabilities, dims=2))
        U .= U * diagm(p) * U'
    end
end

"""
    ProductMeasure(dim::Int)

A type representing a measure on the density states.
It is a product Haar measure on the unitary group and a uniform (Lebesgue) measure on the simplex.
"""
struct ProductMeasure  <: QuantumMeasurementRandom
    dim::Int
end

function Base.rand(rng::AbstractRNG, type::ProductMeasure, n_samples::Integer)
    dim = type.dim
    ps = rand(rng, Simplex(dim), n_samples)
    Us = rand(rng, HaarUnitary(dim), n_samples)
    combine!(Us, ps)
    Us
end

"""
    GinibreEnsamble(dim::Int)

A type representing a Ginibre ensamble of complex matrices of dimension `dim`.
"""
struct GinibreEnsamble  <: QuantumMeasurementRandom
    dim::Int
end

function Base.rand(rng::AbstractRNG, type::GinibreEnsamble, n_samples::Integer)
    ρs = randn(rng, ComplexF32, type.dim, type.dim, n_samples)

    for ρ ∈ eachslice(ρs, dims=3)
        ρ .= ρ * ρ'
        ρ ./= tr(ρ)
    end

    ρs
end

function Base.rand(type::QuantumMeasurementRandom, n_samples::Integer)
    rand(Random.default_rng(), type, n_samples)
end

function Base.rand(rng, type::QuantumMeasurementRandom)
    s = rand(rng, type, 1)
    dropdims(s, dims=ndims(s))
end

function Base.rand(type::QuantumMeasurementRandom)
    s = rand(type, 1)
    dropdims(s, dims=ndims(s))
end