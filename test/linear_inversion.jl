measurements = Dict(
    "Matrix" => assemble_measurement_matrix(polarization_state(s) / √3 for s ∈ (:H, :V, :D, :A, :R, :L)),
    "ProportionalMeasurement Complete" => ProportionalMeasurement(polarization_state(s) for s ∈ (:H, :V, :D, :A, :R, :L)),
    "ProportionalMeasurement Incomplete" => ProportionalMeasurement(polarization_state(s) for s ∈ (:H, :V, :D, :R)),
)

ρs = rand(GinibreEnsamble(2), 10)

for (name, μ) ∈ measurements
    @testset "Exact Linear Inversion ($name)" begin
        for ρ ∈ eachslice(ρs, dims=3)
            θ = traceless_vectorization(ρ)
            freqs = get_probabilities(μ, θ)
            ρ_pred = estimate_state(freqs, μ, LinearInversion())[1]
            @test ρ ≈ ρ_pred
        end
    end

    @testset "Noisy Linear Inversion ($name)" begin
        for ρ ∈ eachslice(ρs, dims=3)
            θ = traceless_vectorization(ρ)
            freqs = get_noisy_probabilities(μ, θ, 0.1)
            ρ_pred = estimate_state(freqs, μ, LinearInversion())[1]
            @test fidelity(ρ, ρ_pred) > 0.99
        end
    end
end