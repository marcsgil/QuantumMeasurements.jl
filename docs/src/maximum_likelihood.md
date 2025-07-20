# Maximum Likelihood Estimation

For situations with limited observations (e.g., photocount regime), the maximum likelihood estimator (MLE) is more appropriate than linear inversion. While linear inversion relies on frequency estimates that become unreliable with few observations, MLE works directly with the count data.

## Theory and Motivation

The maximum likelihood estimator seeks the density operator that maximizes the probability of observing the experimental data:

```math
\hat{\rho}_{ML} = \arg\max_{\rho} \mathcal{L}(N_1, \ldots, N_M | \rho)
```

where the likelihood function for count data is:

```math
\mathcal{L}(N_1, \ldots, N_M | \rho) = \prod_{m=1}^M \left[\mathrm{Tr} \left(\Pi_m \rho \right)\right]^{N_m}
```

Here, $N_m$ represents the number of times outcome $m$ was observed in the experiment, and the total number of observations is $N = \sum_m N_m$.

## Key Advantages

**Statistical Optimality**: MLE is asymptotically efficient, meaning it achieves the Cramér-Rao bound for large sample sizes, providing the best possible precision given the available data.

**Finite Sample Performance**: Unlike linear inversion, which can produce unphysical results (non-positive matrices) when probabilities are poorly estimated from few counts, MLE directly optimizes over valid density operators, ensuring physical outputs.

**Photocount Compatibility**: MLE naturally handles scenarios with discrete photon counts, making it ideal for single-photon experiments and quantum optics applications where the Poisson statistics of detection events matter.

## Implementation Details

The package implements MLE using the **Accelerated Projected Gradient (APG)** algorithm with adaptive restart, following the approach of [PhysRevA.95.062336](@cite). The optimization procedure:

1. **Parameterization**: The density operator is parameterized using the Bloch vector representation to ensure the trace condition $\mathrm{Tr}(\rho) = 1$.

2. **Positivity Constraint**: During optimization, each step checks if the resulting state remains positive semi-definite. If not, the step is rejected and a smaller step size is used.

3. **Convergence**: The algorithm continues until the gradient norm falls below a specified tolerance or the maximum number of iterations is reached.

## Usage

```julia
method = MaximumLikelihood()

ρ_est, θ_est = estimate_state(outcomes, μ, method;
                              x₀=zeros(d²-1),     # Initial Bloch vector
                              t=0.4,              # Initial step size  
                              β=0.8)              # Backtracking factor
```

**Parameters**:
- `x₀`: Initial guess for the Bloch vector (default: maximally mixed state)
- `t`: Initial step size for backtracking line search
- `β`: Reduction factor for backtracking line search

## When to Use MLE

MLE is particularly recommended for:
- **Low photon count experiments** (N < 1000 per measurement setting)
- **Single-photon tomography** where Poisson noise dominates
- **High-precision applications** where statistical optimality is crucial
- **Scenarios where physical constraints must be strictly enforced**

For large sample sizes, linear inversion becomes computationally more efficient and produces similar results. For practical guidance on when to use each method, see [Choosing Estimation Methods](@ref).

## Relationship to Other Methods

MLE provides a middle ground between [Linear Inversion](@ref linear_inversion) and full [Bayesian Inference](@ref):
- **vs Linear Inversion**: More robust for small samples, always produces physical states
- **vs Bayesian Inference**: Faster computation, but no uncertainty quantification

For quantum tomography applications, MLE represents the gold standard for point estimation when computational resources are limited but statistical rigor is required.

## See Also

- [Mathematical Foundations](@ref): Basic theory and linear inversion methods
- [Bayesian Inference](@ref): Full posterior estimation with uncertainty quantification
- [Choosing Estimation Methods](@ref): Practical guide for method selection