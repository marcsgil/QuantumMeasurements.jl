# [Proportional Measurements](@id prop_meas)

Proportional measurements arise when the experimental setup cannot implement a complete POVM, or when the POVM is not properly normalized. Instead, we can only specify a set of positive operators $\{\Pi_m\}$ such that the probability of measuring outcome $m$ is proportional to $\mathrm{Tr}(\rho \Pi_m)$, but the normalization may be state-dependent or unknown. This scenario occurs in several important contexts:

1. **Incomplete measurement settings**: When only a subset of a complete POVM can be measured
2. **Obstructed detection**: When physical obstructions block part of the measurement
3. **Experimental constraints**: When measurement apparatus limitations prevent complete projective measurements

## Mathematical Framework

### The Proportionality Problem

Consider a set of positive operators $\{\Pi_m\}$ where $\sum_m \Pi_m \neq \mathbb{1}$. The measurement yields outcomes with probabilities:

```math
p_m = \frac{\mathrm{Tr}(\rho \Pi_m)}{\mathrm{Tr}(\rho G)}
```

where $G = \sum_m \Pi_m$ is the total measurement operator. The denominator depends on the unknown state $\rho$, creating a normalization problem.

### Transformation to Valid POVM

The `ProportionalMeasurement` type solves this by transforming the problem into standard tomography. The key insight is to:

1. **Compute the Kraus operator**: Find $A$ such that $A^\dagger A = G$ using Cholesky decomposition
2. **Transform measurement operators**: Define $\tilde{\Pi}_m = A^{-\dagger} \Pi_m A^{-1}$
3. **Transform the state**: Work with the transformed state $\tilde{\rho} = A \rho A^\dagger / \mathrm{Tr}(\rho G)$

This transformation ensures $\sum_m \tilde{\Pi}_m = \mathbb{1}$, creating a valid POVM for the transformed state.

### Implementation Details

The [`ProportionalMeasurement`](@ref) constructor automatically handles this transformation:

```julia
# Create from operators that don't form a complete POVM
operators = [π₁, π₂, ..., πₘ]  # where sum(operators) ≠ I
μ = ProportionalMeasurement(operators)

# Use with any estimation method
ρ_est = estimate_state(outcomes, μ, method)[1]
```

The implementation:
- Computes `G = sum(operators)` and its Cholesky decomposition
- Stores the Kraus operator and its inverse for efficient transformations
- Automatically applies the inverse transformation during state reconstruction

## Applications and Examples

### Reduced Measurement Settings

Some experiments use minimal measurement settings for practical reasons. For example, the twin-photon example uses only 16 polarization measurement combinations instead of the complete 36-setting tomography.

**Advantage**: Faster data collection with fewer experimental configurations
**Trade-off**: Requires proportional measurement analysis instead of direct POVM inversion

## When to Use Proportional Measurements

**Use proportional measurements when**:
- Complete POVM implementation is experimentally challenging
- Physical obstructions are unavoidable
- Measurement time must be minimized

**Limitations**:
- Slightly reduced precision compared to complete measurements
- Computational overhead for Kraus transformations
- Requires careful handling of numerical stability in matrix inversions

The `ProportionalMeasurement` framework provides a solution for quantum tomography in realistic experimental conditions where ideal complete measurements are not feasible.