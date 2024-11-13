abstract type AbstractLinearInversion end

"""
    LinearInversion <: AbstractLinearInversion

A type representing a linear inversion method for estimating the state of a quantum system.
The equations are solved by a QR decomposition (`\` operator).
"""
struct LinearInversion <: AbstractLinearInversion end

"""
    estimate_state(outcomes, μ, method; kwargs...)

Estimate the state of a quantum system given the `outcomes` of a set of measurements represented by `μ` using the specified `method`.

If `outcomes` is an array of integers, it is interpreted as the number of times each measurement outcome was observed.
If `outcomes` is an array of floats, it is interpreted as the frequency of each measurement outcome, so the sum of the elements must be 1.

Outputs a tuple `(ρ, θ)` where `ρ` is the estimated density matrix and `θ` is the corresponding Bloch vector.
Additionally, the [`BayesianInference`](@ref) method also returns a covariance matrix `Σ` 
representing the uncertainty in the estimated Bloch vector.

For `AbstractLinearInversion` ([`LinearInversion`](@ref), [`PreAllocatedLinearInversion`](@ref) and [`NormalEquations`](@ref)) methods, 
there are no additional keyword arguments.

For the [`MaximumLikelihood`](@ref) method, the following keyword arguments are available:

- `x₀=zeros(get_measurement_type(μ), get_dim(μ)^2 - 1)`: Initial guess for the Bloch vector. Default is a maximally mixed state.
- `t=0.4`: Initial step size for the backtracking line search.
- `β=0.8`: Factor for the backtracking line search.
- `max_iter=10^3`: Maximum number of iterations.
- `tol=1e-10`: Tolerance for the convergence criterion.

For the [`BayesianInference`](@ref) method, the following keyword arguments are available:

- `verbose=false`: whether to print information about the run.
- `σ=get_measurement_type(μ)(1e-2)`: the initial standard deviation of the proposal distribution.
- `log_prior=θ -> zero(get_measurement_type(μ))`: the log-prior function.
- `θ₀=zeros(get_measurement_type(μ), get_dim(μ)^2 - 1)`: the initial state of the chain.
- `nsamples=10^4`: the number of samples to take.
- `nwarm=10^3`: the number of warm-up samples to take.
- `chain=nothing`: if not `nothing`, store the chain in this matrix.
"""
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

"""
    PreAllocatedLinearInversion(μ) <: AbstractLinearInversion

Create a `PreAllocatedLinearInversion` instance for the given measurement `μ`.

This method precomputes the pseudo-inverse of the corresponding measurement matrix, so that subsequent calls to `estimate_state` are faster.
"""
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


"""
    NormalEquations(μ) <: AbstractLinearInversion

Create a `NormalEquations` instance for the given measurement `μ`.
This method solves the normal equations for the given measurement matrix. This can be faster for systems with a large number of measurements.
"""
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