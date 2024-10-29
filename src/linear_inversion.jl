abstract type AbstractLinearInversion end

struct LinearInversion <: AbstractLinearInversion end

function estimate_state(freqs::AbstractVector{<:AbstractFloat}, μ, ::LinearInversion)
    q = freqs .- get_trace_part(μ) ./ √get_dim(μ)

    θ = get_traceless_part(μ) \ q |> Array
    ρ = density_matrix_reconstruction(θ)
    project2density!(ρ)
    traceless_vectorization!(θ, ρ)

    post_estimation_routine!(ρ, θ, μ)

    ρ, θ
end

function estimate_state(freqs::AbstractArray{<:AbstractFloat}, μ, mthd::AbstractLinearInversion)
    estimate_state(vec(freqs), μ, mthd)
end

function estimate_state(freqs::AbstractArray{<:Integer}, μ, mthd::AbstractLinearInversion)
    estimate_state(normalize(freqs, 1), μ, mthd)
end

struct PreAllocatedLinearInversion{T1<:AbstractMatrix,T2<:AbstractVector} <: AbstractLinearInversion
    pseudo_inv::T1
    θ_correction::T2
end

function PreAllocatedLinearInversion(μ)
    pseudo_inv = pinv(get_traceless_part(μ))
    θ_correction = pseudo_inv * get_trace_part(μ)
    PreAllocatedLinearInversion{typeof(pseudo_inv),typeof(θ_correction)}(pseudo_inv, θ_correction)
end

function estimate_state(freqs::AbstractVector{<:AbstractFloat}, μ, mthd::PreAllocatedLinearInversion)
    θs = similar(mthd.θ_correction)
    T = eltype(θs)
    copy!(θs, mthd.θ_correction)
    mul!(θs, mthd.pseudo_inv, freqs, one(T), -one(T) / √get_dim(μ))
    θs = Array(θs)

    ρ = density_matrix_reconstruction(θs)
    project2density!(ρ)
    post_estimation_routine!(ρ, θs, μ)

    ρ, θs
end

struct NormalEquations{T1<:AbstractMatrix,T2<:AbstractVector} <: AbstractLinearInversion
    TdagT::T1
    Tdagq::T2
end

function NormalEquations(μ)
    T = get_traceless_part(μ)
    TdagT = similar(T, size(T, 2), size(T, 2))
    Tdagq = similar(T, size(T, 2))
    NormalEquations{typeof(TdagT),typeof(Tdagq)}(TdagT, Tdagq)
end

function estimate_state(freqs::AbstractVector{<:AbstractFloat}, μ, method::NormalEquations)
    traceless_part = get_traceless_part(μ)
    mul!(method.TdagT, traceless_part', traceless_part)
    mul!(method.Tdagq, traceless_part', freqs .- get_trace_part(μ) ./ √get_dim(μ))
    θs = method.TdagT \ method.Tdagq |> Array
    ρ = density_matrix_reconstruction(θs)
    project2density!(ρ)
    post_estimation_routine!(ρ, θs, μ)
    ρ, θs
end