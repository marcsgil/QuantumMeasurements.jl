### Measurement interface represantation for AbstractMatrix

"""
    get_traceless_part(μ)

Returns the traceless part of the measurement `μ`.
"""
get_traceless_part(μ::AbstractMatrix) = @view μ[:, begin+1:end]

"""
    get_trace_part(μ)

Returns the trace part of the measurement `μ`.
"""
get_trace_part(μ::AbstractMatrix) = @view μ[:, begin]

"""
    get_dim(μ)

Returns the dimension of the measurement `μ`.
"""
function get_dim(μ::AbstractMatrix)
    dim = √size(μ, 2)
    if !isinteger(dim)
        throw(ArgumentError("The provided AbstractMatrix does not represent a valid measurement. 
    Its second axes has size $(size(μ, 2)), which is not a perfect square."))
    else
        return Int(dim)
    end
end

"""
    get_num_outcomes(μ)

Returns the number of possible outcomes for a single trial of measurement `μ`.
"""
get_num_outcomes(μ::AbstractMatrix) = size(μ, 1)

"""
    get_measurement_type(μ)

Returns the type of the measurement `μ`.
"""
get_measurement_type(μ::AbstractMatrix) = eltype(μ)

"""
    post_estimation_routine!(ρ, θ, μ)

Performs a post-estimation routine on the state `ρ`, `θ` for the measurement `μ`.
Default implementation does nothing.
"""
post_estimation_routine!(ρ, θ, μ) = nothing

"""
    get_probabilities!(dest, μ, θ)

Returns the probabilities of the outcomes of the measurement `μ` given the coefficients `θ`.
The result is stored in `dest`.
"""
function get_probabilities!(dest, μ::AbstractMatrix, θ)
    dim = get_dim(μ)
    trace_part = get_trace_part(μ)
    traceless_part = get_traceless_part(μ)

    T = eltype(dest)

    copy!(dest, trace_part)
    mul!(dest, traceless_part, θ, one(T), convert(T, 1 / √dim))
end

"""
    filter_measurement(μ, J)

Filters the measurement `μ` so that only the outcomes `J` are represented.
"""
filter_measurement(μ::AbstractMatrix, J) = μ[J, :]

"""
    empty_measurement(num_outcomes, dim, T)

Returns an empty measurement matrix with `num_outcomes` outcomes and dimension `dim`.

The type T should be either an `AbstractMatrix` or a `ProportionalMeasurement`(@ref).
"""
function empty_measurement(num_outcomes, dim, ::Type{T}) where {T<:AbstractMatrix}
    T(undef, num_outcomes, dim^2)
end

"""
    update_measurement!(μ, itr)
    update_measurement!(μ::AbstractMatrix, buffer, itr, pars, f!)

Updates the measurement `μ` with an iterator `itr`.
The second sigature is used when the measurement is updated with a buffer and a function `f!`, which is used to update the buffer.
f! should have the signature `f!(buffer, pars)`.
"""
function update_measurement!(μ::AbstractMatrix, itr)
    for (row, Π) ∈ zip(eachrow(μ), itr)
        vectorization!(row, Π)
    end
end

function update_measurement!(μ::AbstractMatrix, buffer, itr, pars, f!)
    for (row, r) ∈ zip(eachrow(μ), itr)
        f!(buffer, r, pars)
        vectorization!(row, buffer)
    end
end

"""
    multithreaded_update_measurement!(μ, buffers, itr, pars, f!)

Same as `update_measurement!`(@ref), but uses multithreading.
"""
function multithreaded_update_measurement!(μ::AbstractMatrix, buffers, itr, pars, f!)
    num_buffers = size(buffers, ndims(buffers))
    chunk_size = cld(length(itr), num_buffers)

    idxs_chunks = Iterators.partition(axes(μ, 1), chunk_size)
    buffers_chunks = eachslice(buffers, dims=ndims(buffers))
    itr_chunks = Iterators.partition(itr, chunk_size)

    tasks = map(zip(idxs_chunks, buffers_chunks, itr_chunks)) do (idxs, buffer, itr)
        @spawn update_measurement!(view(μ, idxs, :), buffer, itr, pars, f!)
    end

    fetch.(tasks)
    nothing
end

### The interface ends here

"""
    get_probabilities(μ, θ)

Returns the probabilities of the outcomes of the measurement `μ` given the coefficients `θ`.
"""
function get_probabilities(μ, θ)
    dest = Vector{get_measurement_type(μ)}(undef, get_num_outcomes(μ))
    get_probabilities!(dest, μ, θ)
end

"""
    assemble_measurement_matrix(itr)

Assembles a measurement matrix from an iterator `itr`.
Each row of the matrix is a vectorization of the corresponding element of `itr`.
`vectorization` is implemented when the element of `itr` is a matrix or a vector, representing a
Positive Operator Valued Measure (POVM) or a Projection Valued Measure (PVM) respectively.
Given `μ = assemble_measurement_matrix(itr)` and a Bloch vector `θ`, the probabilities of the outcomes are calculated by `μ * θ`
"""
function assemble_measurement_matrix(itr)
    hcat((vectorization(Π) for Π ∈ itr)...)'
end