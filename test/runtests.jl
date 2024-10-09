using QuantumMeasurements, Random, Test, LinearAlgebra, Distributions
Random.seed!(123)

function simulate_outcomes(state, measurement, N; atol=1e-4)
    probs = get_probabilities(measurement, traceless_vectorization(state))
    s = sum(probs)

    @assert isapprox(s, 1; atol) """\n The probabilities do not sum to 1, but to $s.
        If you believe this is due to numerical errors, you can try to increase the `atol` parameter.
        """

    normalize!(probs, 1)
    rand(Multinomial(N, probs))
end

include("vec_and_rec.jl")
include("pvmXpovm.jl")
include("linear_inversion.jl")
include("tomography.jl")