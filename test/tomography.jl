methods = Dict(
    "Linear Inversion" => LinearInversion(),
    "Maximum Likelihood" => MaximumLikelihood(),
    "Bayesian Inference" => BayesianInference()
)

Ns = Dict(
    "Linear Inversion" => 10^6,
    "Maximum Likelihood" => 10^4,
    "Bayesian Inference" => 10^4
)

representations = Dict(
    "Projective" => x -> x / √3,
    "Matrix" => x -> x * x' / 3
)

pol_states = Dict(
    "Complete" => [:H, :V, :D, :A, :R, :L],
    "Partial" => [:H, :V, :D, :R]
)

constructors = Dict(
    "Complete" => assemble_measurement_matrix,
    "Partial" => ProportionalMeasurement
)


function test_method(method, N, ρ::AbstractMatrix, measurement)
    outcomes = simulate_outcomes(ρ, measurement, N)
    σ = estimate_state(outcomes, measurement, method)[1]
    @test fidelity(ρ, σ) ≥ 0.99
    @test isposdef(σ)
end

function test_method(method, N, ψ::AbstractVector, measurement)
    outcomes = simulate_outcomes(ψ * ψ', measurement, N)
    σ = estimate_state(outcomes, measurement, method)[1]
    φ = project2pure(σ)
    @test fidelity(ψ, φ) ≥ 0.99
end

for method ∈ keys(methods), representation ∈ keys(representations), constructor ∈ keys(constructors)
    @testset "$method; $representation; $constructor" begin
        for ψ ∈ eachslice(rand(HaarVector(2), 10), dims=2)
            measurement = constructors[constructor](representations[representation](polarization_state(s)) for s in pol_states[constructor])
            test_method(methods[method], Ns[method], ψ, measurement)
        end

        for ρ ∈ eachslice(rand(GinibreEnsamble(2), 10), dims=3)
            measurement = constructors[constructor](representations[representation](polarization_state(s)) for s in pol_states[constructor])
            test_method(methods[method], Ns[method], ρ, measurement)
        end
    end
end