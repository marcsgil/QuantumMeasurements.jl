@setup_workload begin
    symbols1 = [:H, :V, :D, :A, :R, :L]
    symbols2 = [:H, :V, :D, :R]

    @compile_workload begin
        μ1 = assemble_measurement_matrix(polarization_state(s) / √3 for s in symbols1)
        μ2 = ProportionalMeasurement(polarization_state(s) for s in symbols2)

        methods = [LinearInversion(), MaximumLikelihood(), BayesianInference()]
        states = [rand(GinibreEnsamble(2)), rand(HaarVector(2))]

        for method ∈ methods, state ∈ states, μ ∈ [μ1, μ2]
            outcomes = get_probabilities(μ, traceless_vectorization(state))
            σ = estimate_state(outcomes, μ, method)[1]
            fidelity(state, σ)
        end
    end
end