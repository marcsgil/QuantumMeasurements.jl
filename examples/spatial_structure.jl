# # Spatial structure of light

# This example demonstrates how to perform the tomography of the spatial structure of light.
# Although this is not necessarily quantum, the methods provided by 
# QuantumMeasurements can be used to perform this task.
# For this example, we will also use the package 
# [StructuredLight.jl](https://github.com/marcsgil/StructuredLight.jl) for 
# the generation of the structured light patterns, as well as
# [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) for plotting.

using QuantumMeasurements, LinearAlgebra, StructuredLight, CairoMakie, Random
Random.seed!(1234);

# One of the most well-known examples of structured light is the Laguerre-Gauss modes.
# These modes are solutions of the paraxial wave equation in cylindrical coordinates.
# They are characterized by two indices, $p$ and $l$, the first one being the radial index, 
# while the second one is the azimuthal index.
# We can generate these modes using the `lg` function from the StructuredLight package.
rs = LinRange(-3, 3, 256)
ϕ₁ = lg(rs, rs, p=1)
ϕ₂ = lg(rs, rs, l=2)

visualize(abs2.(cat(ϕ₁, ϕ₂, dims=3)))

# We can now create a random superposition of these two modes.
ψs = randn(ComplexF64, 2)
normalize!(ψs)
ψ = ψs[1] * ϕ₁ + ψs[2] * ϕ₂
probs = abs2.(ψ)
normalize!(probs, 1)
visualize(probs)

# What we want to do is recover the coefficients `ψs` from the image.
# We show now that this is equivalent to quantum state tomography.

# ## Theory

# We first notice that the modes are normalized, i.e., $\int |\psi|^2 dA = 1$:
δA = step(rs)^2

isapprox(sum(abs2, ϕ₁) * δA, 1, atol=1e-5) &&
isapprox(sum(abs2, ϕ₂) * δA, 1, atol=1e-5) &&
isapprox(sum(abs2, ψ) * δA, 1, atol=1e-5)

# Therefore, we can regard $|\psi|^2$ as a probability density function.
# The probability of making a detection in region $R$ is given by
# ```math
# p(R) = \int_R |\psi(\mathbf{r})|^2 dA \approx \delta A \times |\psi(\mathbf{r_R})|^2
# ```
# where $\mathbf{r_R}$ is any point in $R$ and we have made the approximation that 
# the region $R$ is small enough that the mode can be considered constant.
# This is a reasonable approximation for pixels in a camera, for example.
# By expanding the mode in the basis $\{\phi_1, \phi_2\}$, we have
# ```math
# p(R) ≈ \sum_{i,j} ψ_i^* ψ_j \phi_i^*(\mathbf{r}_R) \phi_j(\mathbf{r}_R) \delta A = |\braket{\psi|\phi_R}|^2
# ```
# where $\phi_R = (\phi_1^*(\mathbf{r}_R), \ \phi_2^*(\mathbf{r}_R)) \times \sqrt{\delta A}$.
# We see then that the detection probability can be written as the Born rule for a projective measurement in $\phi_R$.

# ## Tomography

# By having understood that this is a quantum state tomography problem, 
# we can proceed to estimate the state.
# We assemble our measurement matrix:
μ = assemble_measurement_matrix([conj.(pair)...] * √δA for pair ∈ zip(ϕ₁, ϕ₂));

# Finally, the measurements are simply the (normalized) intensity.
method = LinearInversion()
ρ = estimate_state(probs, μ, method)[1]

# We can now compare the estimated state with the true state.
fidelity(ρ, ψs) ≈ 1
# We see that the fidelity is close to 1, meaning that we have successfully recovered the state.