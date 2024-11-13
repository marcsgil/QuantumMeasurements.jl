# Quick start

In this tutorial, we will lay down the first steps necessary to perform quantum state tomography. We will learn how to specify what measurement was made, how to represent the experimental results and, finally, how to actually get an estimation of the density operator that produced such experimental results. As an example, we will perform the tomography of the polarization state of a single photon.

As a quick reference, the entire code for this tutorial is shown bellow. We will, in what follows, go through it step by step.

```jldoctest
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

# output

2×2 Matrix{ComplexF64}:
 1.0+0.0im                 0.0+1.4417e-25im
 0.0-1.4417e-25im  5.55112e-17+0.0im
```

## Specifying our measurement

First, we import the package by running `using QuantumMeasurements`. Then, we specify the polarization basis in which the measurements are performed. In this case, we have chosen the horizontal (H), vertical (V), diagonal (D), anti-diagonal (A), right circular (R) and left circular (L) polarization states. This is run by calling
```julia
measurement_basis = [polarization_state(s) / √3 for s ∈ (:H, :V, :D, :A, :R, :L)]
```
The `polarization_state` function returns the polarization state corresponding to the symbol `s`. These are simply two dimensional complex vectors representing the polarization states. We normalize these states by a factor of `√3` to ensure that our set forms a Projective Valued Measure (PVM). For more information on this, check the [Explanation](explanation.md) section of the documentation.

We then assemble a measurement matrix corresponding to this basis by calling
```julia
μ = assemble_measurement_matrix(measurement_basis)
```
This matrix has the property that, given a Bloch vector $\boldsymbol{\theta}$ that represents a state $\rho$, the probability of obtaining the outcome $i$ is given by $p_i = \sum_{j} \mu_{ij} \theta_j$. This is the representation of our measurement that will be used to estimate the state.

## Specifying our measurement results

The experimental outcomes are represented by a vector of integers. The order corresponds to the one defined by the measurement basis. In this case, we had 100 detections for the H projection, 0 detections for the V projection, and 50 detections for each of the remaining four projections. This is represented by the vector
```julia
outcomes = [100, 0, 50, 50, 50, 50]
```

## Performing the tomography