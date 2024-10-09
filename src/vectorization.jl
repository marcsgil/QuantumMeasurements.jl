# Vectorization of traceless part

"""
    check_traceless_compat(θ[, A])

Check if the vector `θ` and the matrix `A` are compatible for traceless vectorization.
If they are, return the dimension of the corresponding space.
"""
function check_traceless_compat(θ, A)
    dim = size(A, 1)
    @assert dim == size(A, 2) "A must be square"
    if length(θ) != dim^2 - 1
        throw(DimensionMismatch("θ has length $(length(θ)), while A has dimension $dim"))
    end
    dim
end

function check_traceless_compat(θ)
    dims = sqrt(length(θ) + 1)
    if !isinteger(dims)
        throw(DimensionMismatch("θ has length $(length(θ)), which is not a perfect square minus 1"))
    else
        dims = Int(dims)
    end
end

"""
    traceless_vectorization!(θ, A)

Vectorize the traceless part of `A` into a vector `θ` of length `dim(A)^2 - 1`.
"""
function traceless_vectorization!(θ, A)
    dim = check_traceless_compat(θ, A)

    T = eltype(θ)
    sqrt2 = convert(T, sqrt(2))

    for n ∈ eachindex(θ)
        if n ≤ dim * (dim - 1) ÷ 2
            I = nth_off_diagonal(n)
            θ[n] = sqrt2 * real(A[I...])
        elseif n ≤ dim * (dim - 1)
            I = nth_off_diagonal(n - dim * (dim - 1) ÷ 2)
            θ[n] = sqrt2 * imag(A[I...])
        else
            j = n - dim * (dim - 1)
            θ[n] = real((sum(j -> A[j, j], 1:j) - j * A[j+1, j+1]) / convert(T, √(j * (j + 1))))
        end
    end
end

struct Projector{T<:AbstractVector}
    ψ::T
end

function Base.getindex(p::Projector, I::Vararg{Int,2})
    p.ψ[I[1]] * conj(p.ψ[I[2]])
end

Base.size(p::Projector) = (length(p.ψ), length(p.ψ))
Base.size(p::Projector, i::Int) = i ≤ 2 ? length(p.ψ) : 1
LinearAlgebra.tr(p::Projector) = sum(abs2, p.ψ)

traceless_vectorization!(θ, ψ::AbstractVector) = traceless_vectorization!(θ, Projector(ψ))

function traceless_vectorization(A)
    θ = Vector{real(eltype(A))}(undef, size(A, 1)^2 - 1)
    traceless_vectorization!(θ, A)
    θ
end

function diag_retrival(n, θz)
    result = zero(eltype(θz))
    for j ∈ eachindex(θz)
        factor = if n ≤ j
            1
        elseif n == j + 1
            -j
        else
            0
        end

        result += θz[j] * factor / √(j^2 + j)
    end
    result
end

function traceless_reconstruction!(A, θ)
    dim = check_traceless_compat(θ, A)

    T = eltype(A)
    inv_sqrt2 = convert(T, inv(sqrt(2)))

    n = 1
    # Off-diagonal elements (real part)
    for k ∈ 1:dim*(dim-1)÷2
        i, j = nth_off_diagonal(k)
        A[i, j] = θ[n] * inv_sqrt2
        A[j, i] = θ[n] * inv_sqrt2
        n += 1
    end

    # Off-diagonal elements (imaginary part)
    for k ∈ 1:dim*(dim-1)÷2
        i, j = nth_off_diagonal(k)
        A[i, j] += im * θ[n] * inv_sqrt2
        A[j, i] -= im * θ[n] * inv_sqrt2
        n += 1
    end

    θz = @view θ[dim*(dim-1)+1:end]
    # Diagonal elements
    for j in 1:dim
        A[j, j] = diag_retrival(j, θz)
    end
end

function traceless_reconstruction(θ)
    dims = check_traceless_compat(θ)

    A = Matrix{complex(eltype(θ))}(undef, dims, dims)
    gell_mann_reconstruction!(A, θ)
    A
end

# Vectorization of density matrix

"""
    density_matrix_reconstruction!(ρ, θ)

Reconstruct a density matrix from the vector `θ` and stores the result in `ρ`.
"""
function density_matrix_reconstruction!(ρ, θ)
    traceless_reconstruction!(ρ, θ)
    dim = size(ρ, 1)
    factor = convert(eltype(ρ), 1 / dim)
    for i ∈ 1:dim
        ρ[i, i] += factor
    end
end

"""
    density_matrix_reconstruction(θ)

Reconstruct a density matrix from the vector `θ`.
"""
function density_matrix_reconstruction(θ)
    dims = check_traceless_compat(θ)
    ρ = Matrix{complex(eltype(θ))}(undef, dims, dims)
    density_matrix_reconstruction!(ρ, θ)
    ρ
end

# Vectorization of hermitian operators

function vectorization!(x, A)
    x[begin] = real(tr(A)) / √size(A, 1)
    traceless_vectorization!((@view x[begin+1:end]), A)
end

vectorization!(x, ψ::AbstractVector) = vectorization!(x, Projector(ψ))

function vectorization(M)
    T = real(eltype(M))
    dim = size(M, 1)
    x = Vector{T}(undef, dim^2)
    vectorization!(x, M)
    x
end

function reconstruction!(A, x)
    traceless_reconstruction!(A, (@view x[begin+1:end]))
    dim = size(A, 1)
    factor = convert(eltype(A), 1 / sqrt(dim))
    for n ∈ 1:dim
        A[n, n] += x[begin] * factor
    end
end

function reconstruction(x)
    dim = Int(sqrt(length(x)))
    M = Matrix{complex(eltype(x))}(undef, dim, dim)
    reconstruction!(M, x)
    M
end