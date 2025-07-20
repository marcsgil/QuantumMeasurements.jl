# Choosing Estimation Methods

This guide helps you select the appropriate quantum state tomography method for your specific experimental conditions and requirements.

## Quick Decision Tree

```
Do you have complete POVM measurements?
├─ No → Use ProportionalMeasurement wrapper and continue below
└─ Yes → Continue below

How many photon counts do you have?
├─ High (>1000 per setting) → Linear methods recommended
│  ├─ Single estimation → LinearInversion
│  ├─ Many estimations (same setup) → PreAllocatedLinearInversion  
│  └─ Highly redundant measurements → NormalEquations
└─ Low (<1000 per setting) → Advanced methods recommended
   ├─ Need uncertainty quantification → BayesianInference
   └─ Want optimal point estimate → MaximumLikelihood
```

## Method Comparison

| Method | Speed | Data Requirements | Output | Best For |
|--------|-------|------------------|---------|----------|
| `LinearInversion` | Fast | High counts | Point estimate | General purpose, high SNR |
| `PreAllocatedLinearInversion` | Very Fast* | High counts | Point estimate | Repeated measurements |
| `NormalEquations` | Fast | High counts | Point estimate | Overcomplete measurements |
| `MaximumLikelihood` | Medium | Any counts | Point estimate | Low counts, optimal accuracy |
| `BayesianInference` | Very Slow | Any counts | Full posterior | Uncertainty quantification |

*After preprocessing overhead

## Detailed Method Selection

### Linear Inversion Methods

**When to use**: High count experiments.

#### LinearInversion
- **Use for**: Single state estimations or changing measurement setups
- **Advantages**: Simple, robust, numerically stable
- **Limitations**: No uncertainty quantification, can produce unphysical states with low counts

```julia
method = LinearInversion()
ρ, θ = estimate_state(outcomes, measurement, method)
```

#### PreAllocatedLinearInversion  
- **Use for**: Repeated tomography with the same measurement setup
- **Advantages**: Significant speedup for batch processing, same accuracy as LinearInversion
- **Limitations**: Memory overhead, setup cost, inflexible to measurement changes

```julia
method = PreAllocatedLinearInversion()
# Use same method instance for multiple estimations
ρ₁, θ₁ = estimate_state(outcomes₁, measurement, method)
ρ₂, θ₂ = estimate_state(outcomes₂, measurement, method)  # Much faster
```

#### NormalEquations
- **Use for**: Measurements with many more POVM elements than parameters (M ≫ d²-1)
- **Advantages**: Computational efficiency for overcomplete measurements
- **Limitations**: Less numerically stable, no benefit for moderately sized problems

```julia
method = NormalEquations()
ρ, θ = estimate_state(outcomes, measurement, method)
```

### Advanced Methods

**When to use**: Low photon count experiments, single-photon tomography, or when statistical optimality/uncertainty is crucial.

#### MaximumLikelihood
- **Use for**: Photocount regime, optimal point estimates
- **Advantages**: Handles Poisson statistics correctly, always produces physical states, asymptotically optimal
- **Limitations**: Slower than linear methods, no uncertainty quantification

```julia
method = MaximumLikelihood()
ρ, θ = estimate_state(outcomes, measurement, method;
                      x₀=initial_guess,  # optional
                      t=0.4,             # step size
                      β=0.8)             # backtracking factor
```

#### BayesianInference  
- **Use for**: When uncertainty quantification is essential
- **Advantages**: Full posterior distribution, natural error bars, handles any sample size
- **Limitations**: Very slow (orders of magnitude), requires careful prior selection

```julia
method = BayesianInference()
ρ, θ, Σ = estimate_state(outcomes, measurement, method;
                         nsamples=10000,    # posterior samples
                         nwarm=1000,        # burn-in period
                         σ=1e-2)            # step size
# Σ contains covariance matrix for error bars
```

## Measurement Types

### Standard Measurements
Use when your POVM elements sum to identity: $\sum_m \Pi_m = \mathbb{1}$

```julia
# Create measurement matrix directly
measurement = assemble_measurement_matrix(povm_elements)
```

### Proportional Measurements  
Use when measurements are incomplete or obstructed: $\sum_m \Pi_m \neq \mathbb{1}$

```julia
# Wrapper handles normalization issues automatically
measurement = ProportionalMeasurement(povm_elements)
```

**Common scenarios**:
- Spatial measurements with beam obstructions
- Reduced measurement settings for efficiency
- Incomplete POVM implementations due to experimental constraints

## Practical Guidelines

### High-Throughput Tomography
For real-time or high-volume tomography:
1. Use `PreAllocatedLinearInversion` for speed
2. Ensure sufficient counts (>1000 per setting)
3. Validate occasionally with `MaximumLikelihood`

### Low counts: 
1. Always use `MaximumLikelihood` or `BayesianInference`
2. Consider prior knowledge in Bayesian approach

### Benchmark/Validation Studies
For comparing methods or validating algorithms:
1. Generate test data with [`Random States`](@ref random_states)
2. Compare computational efficiency vs accuracy trade-offs

### Debugging Poor Results
If reconstruction fidelity is low:
1. Check measurement informationally complete (matrix rank)
2. Verify sufficient counts for chosen method
3. Consider measurement calibration issues
4. Try different estimation methods for comparison

Choose methods based on your experimental constraints, computational resources, and accuracy requirements. When in doubt, `LinearInversion` provides a robust starting point for high-count data, while `MaximumLikelihood` handles low-count scenarios effectively.