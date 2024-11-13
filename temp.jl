using QuantumMeasurements

# We specify the polarization basis in which the measurements are performed
# For the reason behind the factor of √3, check the Explanation section of the docs.
measurement_basis = [polarization_state(s) / √3 for s ∈ (:H, :V, :D, :A, :R, :L)]

# We assemble a measurement matrix corresponding to this basis
μ = assemble_measurement_matrix(measurement_basis)

# The experimental outcomes are represented by a vector of integers
# The order corresponds to the one defined by the measurement basis
# In this case, we had 100 detections for the H projection, 0 detections for the V projection,...
outcomes = [100, 0, 50, 50, 50, 50]

# We choose a tomography method
# This is the simplest one
method = LinearInversion()

# Finally, we call `estimate_state` to get an estimate of the state
ρ = estimate_state(outcomes, μ, method)[1]