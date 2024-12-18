struct ProportionalMeasurement{T<:AbstractFloat,TM<:AbstractMatrix{T}}
    measurement_matrix::TM
    measurement_sum::Matrix{Complex{T}}
    kraus_operator::UpperTriangular{Complex{T},Matrix{Complex{T}}}
    inv_kraus_operator::UpperTriangular{Complex{T},Matrix{Complex{T}}}
end

get_traceless_part(μ::ProportionalMeasurement) = get_traceless_part(μ.measurement_matrix)
get_trace_part(μ::ProportionalMeasurement) = get_trace_part(μ.measurement_matrix)
get_dim(μ::ProportionalMeasurement) = get_dim(μ.measurement_matrix)
get_num_outcomes(μ::ProportionalMeasurement) = get_num_outcomes(μ.measurement_matrix)
get_measurement_type(μ::ProportionalMeasurement) = get_measurement_type(μ.measurement_matrix)

function kraus_transformation!(ρ, A)
    rmul!(ρ, A')
    lmul!(A, ρ)
end

function kraus_transformation(ρ, A)
    σ = copy(ρ)
    kraus_transformation!(σ, A)
    σ
end

function kraus_transformation(ψ::AbstractVector, A)
    A * ψ
end

function post_measurement_state!(ρ, A)
    kraus_transformation!(ρ, A)
    ρ ./= tr(ρ)
end

function post_estimation_routine!(ρ, θ, μ::ProportionalMeasurement)
    post_measurement_state!(ρ, μ.inv_kraus_operator)
    project2density!(ρ)
    traceless_vectorization!(θ, ρ)
end

filter_measurement(μ::ProportionalMeasurement, J) = μ.measurement_matrix[J, :]

_transform(ρ::AbstractMatrix) = ρ
_transform(ψ::AbstractVector) = ψ * ψ'

"""
    ProportionalMeasurement(itr)

Constructs a proportional measurement from the iterator `itr`.

This measurement is proportional in the sense that, given Given `μ = ProportionalMeasurement(itr)` 
and a Bloch vector `θ`, the probabilities of the outcomes are proportional to `μ.measurement_matrix * θ`
"""
function ProportionalMeasurement(itr)
    measurement_sum = sum(_transform, itr)
    kraus_operator = cholesky(measurement_sum).U
    inv_kraus_operator = inv(kraus_operator)

    transformed_itr = Iterators.map(Π -> kraus_transformation(Π, inv_kraus_operator'), itr)
    measurement_matrix = assemble_measurement_matrix(transformed_itr)

    ProportionalMeasurement(measurement_matrix, measurement_sum, kraus_operator, inv_kraus_operator)
end

function empty_measurement(num_outcomes, dim, ::Type{ProportionalMeasurement{T,TM}}) where {T,TM}
    measurement_matrix = Matrix{T}(undef, num_outcomes, dim^2)
    measurement_sum = Matrix{complex(T)}(undef, dim, dim)
    kraus_operator = UpperTriangular(measurement_sum)
    inv_kraus_operator = UpperTriangular(measurement_sum)
    ProportionalMeasurement(measurement_matrix, measurement_sum, kraus_operator, inv_kraus_operator)
end

function update_measurement!(μ::ProportionalMeasurement, itr)
    μ.measurement_sum .= sum(_transform, itr)
    μ.kraus_operator .= cholesky(μ.measurement_sum).U
    μ.inv_kraus_operator .= inv(μ.kraus_operator)

    transformed_itr = Iterators.map(Π -> kraus_transformation(Π, μ.inv_kraus_operator'), itr)
    update_measurement!(μ.measurement_matrix, transformed_itr)
end

function update_measurement!(μ::ProportionalMeasurement, buffer, itr, pars, f!)
    μ.measurement_sum .= sum(r -> f!(buffer, r, pars) |> _transform, itr)
    μ.kraus_operator .= cholesky(μ.measurement_sum).U
    μ.inv_kraus_operator .= inv(μ.kraus_operator)
    function transformed_f!(buffer, r, pars)
        f!(buffer, r, pars)
        kraus_transformation(buffer, μ.inv_kraus_operator')
    end
    update_measurement!(μ.measurement_matrix, buffer, itr, pars, transformed_f!)
end

function get_probabilities!(dest, μ::ProportionalMeasurement, θ)
    ρ = density_matrix_reconstruction(θ)
    post_measurement_state!(ρ, μ.kraus_operator)
    traceless_vectorization!(θ, ρ)
    get_probabilities!(dest, μ.measurement_matrix, θ)
end