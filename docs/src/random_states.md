# [Random State Generation](@id random_states)

```@setup setup
using QuantumMeasurements, LinearAlgebra
```

QuantumMeasurements.jl provides utilities for generating random quantum states and related mathematical objects. These are useful for testing tomography algorithms, benchmarking estimation methods, and simulating quantum systems.

## Available Random Distributions

The package implements several important random distributions commonly used in quantum information:

### Haar Random Unitaries

[`HaarUnitary`](@ref) generates random unitary matrices distributed according to the Haar measure, the unique unitarily invariant measure on the unitary group.

```@example setup
# Generate a single 3×3 Haar random unitary
U = rand(HaarUnitary(3))

# Generate multiple samples
U_batch = rand(HaarUnitary(3), 4)  # 3×3×4 array
```

**Implementation**: Uses QR decomposition of complex Gaussian matrices with proper phase normalization to ensure uniform distribution on the unitary group.

### Haar Random State Vectors

[`HaarVector`](@ref) generates random pure quantum states (normalized complex vectors) distributed according to the induced Haar measure.

```@example setup
# Generate a random pure state in dimension 4
ψ = rand(HaarVector(4))

# Multiple pure states
ψ_batch = rand(HaarVector(4), 100)  # 4×100 array
```

This corresponds to the first column of a Haar random unitary matrix.

### Simplex Sampling

[`Simplex`](@ref) generates random probability distributions (positive vectors that sum to 1) uniform on the probability simplex.

```@example setup
# Generate random probabilities for 5 outcomes
p = rand(Simplex(5))
sum(p) ≈ 1  # true

# Batch generation
p_batch = rand(Simplex(5), 50)  # 5×50 array
```

**Algorithm**: Uses the stick-breaking construction for uniform sampling on the simplex.

## Random Density Matrices

The package provides two methods for generating random mixed quantum states:

### Product Measure

[`ProductMeasure`](@ref) generates random density matrices by combining Haar random unitaries with uniform random eigenvalue distributions:

```@example setup
# Generate a random density matrix
ρ = rand(ProductMeasure(3))

# Properties
tr(ρ) ≈ 1         # normalized
all(real.(eigvals(ρ)) .≥ 0)  # positive semidefinite
```

This measure factorizes as the product of the Haar measure on unitaries and the uniform (Lebesgue) measure on the probability simplex. It produces states with typically well-separated eigenvalues.

### Ginibre Ensemble

[`GinibreEnsamble`](@ref) generates random density matrices from the Ginibre ensemble - complex Gaussian random matrices normalized to unit trace:

```@example setup
# Generate using Ginibre ensemble  
ρ = rand(GinibreEnsamble(3))
```

**Construction**: Creates a complex Gaussian matrix G, then forms ρ = GG†/Tr(GG†).

## Implementation Notes

- All random generators use `ComplexF32` for memory efficiency while maintaining sufficient precision
- The implementations are optimized for batch generation when multiple samples are needed
- Random number generation respects @example's standard `AbstractRNG` interface for reproducibility

These random state generators provide the foundation for comprehensive testing and validation of quantum tomography algorithms under realistic conditions.