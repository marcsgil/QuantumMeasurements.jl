# # Spatial structure of light with obstructions

using QuantumMeasurements, StructuredLight, CairoMakie, Random
Random.seed!(1234);

rs = LinRange(-3, 3, 256)
ϕ₁ = lg(rs, rs, p=1)
ϕ₂ = lg(rs, rs, l=2)

visualize(abs2.(cat(ϕ₁, ϕ₂, dims=3)))

ψs = randn(ComplexF64, 2)
ψs /= √sum(abs2, ψs)
ψ = ψs[1] * ϕ₁ + ψs[2] * ϕ₂
visualize(abs2.(ψ))
##

# We now introduce an obstruction in the spatial structure.
# This can be done by multiplying the mode by a mask.

mask = [x > 0 for x ∈ rs, y ∈ rs]
ψ_obstructed = ψ .* mask
visualize(abs2.(ψ_obstructed))

##
μ = assemble_measurement_matrix([conj.(pair)...] * √δA for pair ∈ zip(ϕ₁, ϕ₂));

method = LinearInversion()
ρ = estimate_state(abs2.(ψ) * δA, μ, method)[1]

fidelity(ρ, ψs)