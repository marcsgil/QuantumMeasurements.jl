"""
    nth_off_diagonal(n)

Calculate the row and column of the `n`-th off-diagonal element of a square matrix.
"""
function nth_off_diagonal(n)
    column = ceil(Int, (-1 + √(1 + 8n)) / 2)
    row = n - column * (column - 1) ÷ 2
    row, column + 1
end

"""
    GellMannMatrices(dim, ::Type{T}=ComplexF32) where {T<:Complex}

Generate the generalized Gell-Mann matrices of dimension `dim` with elements of type `T`.
"""
struct GellMannMatrices{T<:Complex}
    dim::Int
end

function GellMannMatrices(dim, ::Type{T}=ComplexF32) where {T<:Complex}
    GellMannMatrices{T}(dim)
end

function Base.iterate(iter::GellMannMatrices{T}, state=1) where {T}
    dim = iter.dim
    state == dim^2 && return nothing

    result = zeros(T, (dim, dim))

    if state ≤ dim * (dim - 1) ÷ 2
        i, j = nth_off_diagonal(state)
        result[i, j] = 1 / √2
        result[j, i] = 1 / √2
    elseif state ≤ dim * (dim - 1)
        i, j = nth_off_diagonal(state - dim * (dim - 1) ÷ 2)
        result[i, j] = im / √2
        result[j, i] = -im / √2
    else
        j = state - dim * (dim - 1)
        factor = 1 / √(j * (j + 1))
        for k ∈ 1:j
            result[k, k] = factor
        end
        result[j+1, j+1] = -j * factor
    end

    result, state + 1
end

Base.IteratorSize(::GellMannMatrices) = Base.HasLength()
Base.IteratorEltype(::GellMannMatrices) = Base.HasEltype()
Base.eltype(::Type{GellMannMatrices{T}}) where {T} = Matrix{T}
Base.length(iter::GellMannMatrices) = iter.dim^2 - 1