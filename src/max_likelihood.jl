struct MaximumLikelihood end

function log_likelihood!(∇ℓπ, buffer1, buffer2, frequencies, μ, y, x)
    get_probabilities!(buffer1, μ, y)

    if any(x -> x < 0, buffer1)
        y .= x
        get_probabilities!(buffer1, μ, y)
    end

    broadcast!(/, buffer2, frequencies, buffer1)
    mul!(∇ℓπ, get_traceless_part(μ)', buffer2)
    broadcast!(log, buffer1, buffer1)
    frequencies ⋅ buffer1
end


function log_likelihood!(buffer, frequencies, μ, x)
    get_probabilities!(buffer, μ, x)
    broadcast!(log, buffer, buffer)
    frequencies ⋅ buffer
end

function update_x!(x, y, ρ, t, ∇ℓπ)
    @. x = y + t * ∇ℓπ
    density_matrix_reconstruction!(ρ, x)
    project2density!(ρ)
    traceless_vectorization!(x, ρ)
end

function gradient_ascent!(x, x_prev, y, buffer1, buffer2, ∇ℓπ, ρ, δ, δ_hat, frequencies, μ, t, β, max_iter, tol)
    θ = 1

    for i in 1:max_iter
        ℓ = log_likelihood!(∇ℓπ, buffer1, buffer2, frequencies, μ, y, x)

        update_x!(x, y, ρ, t, ∇ℓπ)
        @. δ = x - y

        # Backtracking line search
        while log_likelihood!(buffer1, frequencies, μ, x) ≤ ℓ + real(∇ℓπ ⋅ δ - (δ ⋅ δ) / (2t))
            t *= β
            update_x!(x, y, ρ, t, ∇ℓπ)
            @. δ = x - y
            iszero(δ) && return nothing
        end

        @. δ_hat = x - x_prev

        if sum(abs2, δ) < tol
            break
        end

        if δ ⋅ δ_hat < 0
            θ = 1
            x .= x_prev
            y .= x_prev
        else
            new_θ = (1 + sqrt(1 + 4 * θ^2)) / 2
            @. y = x + (θ - 1) * δ_hat / new_θ

            x_prev .= x
            θ = new_θ
        end
    end
end

function estimate_state(outcomes, μ, ::MaximumLikelihood;
    x₀=zeros(get_measurement_type(μ), get_dim(μ)^2 - 1),
    t=0.4,
    β=0.8,
    max_iter=10^3,
    tol=1e-10)

    I = findall(!iszero, vec(outcomes))
    frequencies = normalize(outcomes[I], 1)
    fitered_μ = filter_measurement(μ, I)
    buffer1 = similar(x₀, length(I))
    buffer2 = similar(x₀, length(I))
    ∇ℓπ = similar(x₀)
    x = copy(x₀)
    x_prev = copy(x)
    y = copy(x)
    ρ = density_matrix_reconstruction(x)
    δ = similar(x)
    δ_hat = similar(x)

    gradient_ascent!(x, x_prev, y, buffer1, buffer2, ∇ℓπ, ρ, δ, δ_hat, frequencies, fitered_μ, t, β, max_iter, tol)

    density_matrix_reconstruction!(ρ, y)

    post_estimation_routine!(ρ, y, μ)

    ρ, y
end