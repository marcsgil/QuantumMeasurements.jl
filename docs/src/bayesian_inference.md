# Bayesian Inference

Bayesian inference addresses limitations of both linear inversion and maximum likelihood estimation. While linear inversion can produce unphysical states and MLE yields point estimates without uncertainty quantification, Bayesian methods provide posterior distributions that incorporate uncertainty and prior knowledge.

## Theory and Motivation

Bayesian Mean Estimation (BME) computes a posterior distribution over all possible states and reports the posterior mean as the estimate [Blume-Kohout_2010](@cite):

```math
\hat{\rho}_{\text{BME}} = \int \rho \, \pi(\rho | \mathcal{D}) \, d\rho
```

where $\pi(\rho | \mathcal{D})$ is the posterior distribution given data $\mathcal{D}$. According to Bayes' rule:

```math
\pi(\rho | \mathcal{D}) = \frac{\mathcal{L}(\mathcal{D} | \rho) \pi_0(\rho)}{\mathcal{Z}}
```

Here, $\mathcal{L}(\mathcal{D} | \rho)$ is the likelihood function, $\pi_0(\rho)$ is the prior distribution, and $\mathcal{Z}$ is the normalization constant.

## The Prior Distribution

The choice of prior $\pi_0(\rho)$ is both a strength and limitation of Bayesian inference. The prior encodes assumptions about the quantum state before observing data, which can be valuable when genuine prior knowledge exists (e.g., from physical constraints or previous experiments).

**Default Choice**: The package uses a uniform prior in the Bloch vector representation by default (`log_prior=θ -> zero(...)`), which corresponds to a flat distribution over the traceless part of the density matrix. This choice is uninformative and lets the data dominate the posterior.

**Arbitrariness Concern**: The prior choice is inherently subjective and can influence results, particularly with limited data. Different reasonable priors may yield different estimates. However, as the amount of data increases, the prior's influence diminishes and the posterior concentrates around the maximum likelihood estimate.

**Robustness**: For informationally complete measurements with sufficient data, Bayesian estimates become relatively insensitive to prior choice. The likelihood function dominates the posterior, making the method practically robust despite the theoretical arbitrariness.

## Key Properties

**Full-Rank States**: BME never produces zero eigenvalues. The posterior mean incorporates uncertainty, ensuring all eigenvalues remain positive for reasonable priors.

**Uncertainty Quantification**: The posterior covariance matrix provides error bars. Each parameter's uncertainty can be computed from the posterior distribution.

**Optimal Performance**: Under proper scoring rules, BME minimizes the expected loss for operationally meaningful distance measures between quantum states.

**Small Sample Performance**: BME handles limited datasets by incorporating prior knowledge and regularization through posterior averaging.

## Implementation: MALA Algorithm

The package implements Bayesian inference using the Metropolis-Adjusted Langevin Algorithm (MALA), an MCMC method that samples from the posterior distribution.

MALA proposes new states using the gradient of the log-posterior:

```math
\theta_{\text{prop}} = \theta_{\text{current}} + \frac{\sigma^2}{2} \nabla \log \pi(\theta_{\text{current}} | \mathcal{D}) + \sigma \epsilon
```

where $\epsilon \sim \mathcal{N}(0, I)$ is Gaussian noise and $\sigma$ controls the step size.

Each proposal is accepted with probability:

```math
\alpha = \min\left(1, \frac{\pi(\theta_{\text{prop}} | \mathcal{D})}{\pi(\theta_{\text{current}} | \mathcal{D})} \times \frac{q(\theta_{\text{current}} | \theta_{\text{prop}})}{q(\theta_{\text{prop}} | \theta_{\text{current}})}\right)
```

where $q(\theta' | \theta)$ is the proposal distribution probability density for transitioning from state $\theta$ to $\theta'$. In MALA, this follows a multivariate Gaussian distribution centered at the gradient-adjusted current state. The ratio $\frac{q(\theta_{\text{current}} | \theta_{\text{prop}})}{q(\theta_{\text{prop}} | \theta_{\text{current}})}$ corrects for the asymmetry introduced by the gradient-based proposals, ensuring detailed balance and proper sampling from the posterior.

**Implementation Features**:

1. **Adaptive Step Size**: The algorithm adjusts $\sigma$ to target 57.4% acceptance rate
2. **Constraint Handling**: Proposals leading to non-positive matrices are rejected
3. **Warm-up Phase**: Initial burn-in period for chain convergence

## Computational Cost

**Warning**: Bayesian inference is computationally expensive, requiring orders of magnitude more time than linear inversion or MLE. Typical runs require thousands of MCMC samples, each involving gradient computations and matrix operations. This method should be used when uncertainty quantification justifies the computational overhead.

## Usage

```julia
method = BayesianInference()

ρ_est, θ_est, Σ = estimate_state(outcomes, μ, method; 
                                 nsamples=10000,                    # Posterior samples
                                 nwarm=1000,                        # Warm-up iterations
                                 σ=1e-2,                           # Initial step size
                                 log_prior=θ -> zero(eltype(θ)))    # Prior function
```

**Parameters**:
- `nsamples`: Number of MCMC samples to collect from posterior
- `nwarm`: Number of warm-up iterations for chain equilibration
- `σ`: Initial step size parameter (automatically adapted)
- `log_prior`: Function specifying the log-prior density
- `θ₀`: Initial Bloch vector (default: maximally mixed state)
- `chain`: Optional matrix to store the full MCMC chain

**Output**:
- `ρ_est`: Posterior mean density matrix
- `θ_est`: Posterior mean Bloch vector  
- `Σ`: Posterior covariance matrix (uncertainty quantification)

## When to Use Bayesian Inference

Use Bayesian inference when:
- Uncertainty quantification is required
- Working with limited data where other methods fail
- Parameter correlations are important
- Computational cost is acceptable relative to accuracy requirements

For routine tomography or when computational speed is important, consider linear inversion or MLE instead. See [Choosing Estimation Methods](@ref) for detailed guidance.

## Relationship to Other Methods

Bayesian inference represents the most complete statistical approach:
- **vs [Linear Inversion](@ref linear_inversion)**: Handles any sample size, provides uncertainty quantification
- **vs [Maximum Likelihood Estimation](@ref)**: Full posterior vs point estimate, incorporates prior knowledge
- **Cost**: Significantly higher computational expense

The method is recommended when the additional computational cost is justified by the need for rigorous uncertainty analysis or when working in the challenging regime of very limited experimental data.

## See Also

- [Mathematical Foundations](@ref): Basic theory and linear inversion methods  
- [Maximum Likelihood Estimation](@ref): Optimal point estimation for low counts
- [Choosing Estimation Methods](@ref): Practical guide for method selection
- [Random States](@ref random_states): Generating test data for validation