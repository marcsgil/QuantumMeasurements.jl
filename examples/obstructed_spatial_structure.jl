# # Spatial structure of light with obstructions

# This example is similar to [Spatial structure of light](spatial_structure.md), but now we will
# introduce an obstruction in the spatial structure of light. This will demonstrate how the tomography
# can still be performed even when the spatial structure is not fully accessible.

# We will use the same packages as before:

using QuantumMeasurements, LinearAlgebra, StructuredLight, CairoMakie, Random
Random.seed!(1234);

# We now generate the Laguerre-Gauss modes as before, but we will introduce an obstruction
# in the spatial structure of light. The obstruction will be a horizontal barrier that blocks
# the light in the middle of the image. In the code, this is achieved by introducing coordinates
# `xs` that only cover the right half of the image (x > 0), and `ys` that cover the full height:

xs = LinRange(0, 3, 128)
ys = LinRange(-3, 3, 256)

# We generate the Laguerre-Gauss modes as before, but now using the new coordinates:

ϕ₁ = lg(xs, ys, p=1)
ϕ₂ = lg(xs, ys, l=2)

visualize(abs2.(cat(ϕ₁, ϕ₂, dims=3)))

# We can now create a random superposition of these two modes, as before:

ψs = randn(ComplexF64, 2)
normalize!(ψs)
ψ = ψs[1] * ϕ₁ + ψs[2] * ϕ₂
probs = abs2.(ψ)
normalize!(probs, 1)
visualize(probs)

# Now, we use a `ProportionalMeasurement` to perform the tomography to define our measurement, since,
# now, it is no longer a PVM (Projection Valued Measure) due to the obstruction. This measurement
# type already takes care of the logic of the obstruction, and performs the necessary normalization.
# For that reason, we no longer need to multiply by `δA` (the area of the pixel), as in the previous example.
# Take a look at the [Proportional Measurements](proportional_measurements.md) documentation for more details.

μ = ProportionalMeasurement([conj.(pair)...] for pair ∈ zip(ϕ₁, ϕ₂))

# Finally, everything is the same as before, we can estimate the state using the `LinearInversion` method.

method = LinearInversion()
ρ = estimate_state(probs, μ, method)[1]

# We can now compare the estimated state with the true state.
# The fidelity should be close to 1, meaning that we have successfully recovered the state.

fidelity(ρ, ψs) ≈ 1