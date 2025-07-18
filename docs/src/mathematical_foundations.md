# Mathematical Foundations

This section provides the mathematical foundations of quantum state tomography that underpin the methods implemented in QuantumMeasurements.jl.

## Quantum State Tomography

Let $\mathcal{H}$ be a Hilbert space of dimension $d$, and let $\mathrm{Her}(\mathcal{H}) \subset \mathcal{H}$ be the set of Hermitian operators acting on $\mathcal{H}$. The state of a quantum system is represented by an element $\rho$ of the set of positive semi-definite operators $\mathrm{Pos}(\mathcal{H}) \subset \mathrm{Her}(\mathcal{H})$ such that $\mathrm{Tr} \rho = 1$. Quantum state tomography is the process of reconstructing the density operator $\rho$ from experimental data.

## Observables and Measurements

Observables are represented by elements $A \in \mathrm{Her}(\mathcal{H})$, and their expectation values are given by:

```math
\langle A \rangle = \mathrm{Tr} \left(\rho A\right)
```

A Positive Operator Valued Measure (POVM) is a set of observables $\{\Pi_m\} \subset \mathrm{Pos}(\mathcal{H})$ with the property that:

```math
\sum_m \Pi_m = \mathbb{1}
```

where $\mathbb{1}$ is the identity operator. These operators model the possible outcomes of an experiment: outcome $m$ happens with probability $p_m = \mathrm{Tr} \left(\Pi_m \rho \right)$ according to Born's rule. The POVM conditions ensure that $p_m \geq 0$ and $\sum_m p_m = 1$.

A special case of a POVM is a projective measurement (PVM), where each $\Pi_m$ is a projector, i.e., $\Pi_m = |u_m\rangle\langle u_m|$ for some basis vector $|u_m\rangle \in \mathcal{H}$.

## Parameterization of States

Usually, one chooses a specific basis for the Hilbert space $\mathcal{H}$ to represent the state $\rho$, say $\{|u_1\rangle, \ldots, |u_d\rangle\}$. The density operator can then be expressed in terms of the basis vectors as:

```math
\rho = \sum_{j,k=1}^{d} \rho_{jk} |u_j\rangle\langle u_k|
```

where $\rho_{jk} = \langle u_j | \rho | u_k \rangle$ are the matrix elements of $\rho$ in this basis.

Nonetheless, this matrix representation is redundant because $\rho$ is Hermitian and has unit trace. This means that only $d^2 - 1$ independent real parameters are needed to specify $\rho$. We can choose a traceless basis of Hermitian operators to represent the state uniquely. A convenient choice is the generalized Gell-Mann matrices, which come in three types:

```math
\begin{aligned}
X_{jk} &= \frac{|u_j\rangle\langle u_k| + |u_k\rangle\langle u_j|}{\sqrt{2}} \\
Y_{jk} &= \frac{i(|u_j\rangle\langle u_k| - |u_k\rangle\langle u_j|)}{\sqrt{2}} \\
Z_j &= \frac{1}{\sqrt{j + j^2}} \left( \sum_{r=1}^j |u_r\rangle\langle u_r| - j |u_{j+1}\rangle\langle u_{j+1}| \right)
\end{aligned}
```

where $j = 1, \ldots, d-1$, $k = j+1, \ldots, d$, and $\{|u_1\rangle, \ldots, |u_d\rangle\}$ is our chosen basis set of $\mathcal{H}$. The type [`GellMannMatrices`](@ref) is an iterator over these matrices.

We can specify any state $\rho$ by a list of real coefficients $\boldsymbol{\theta} = (\theta_1,\ldots,\theta_{d^2-1})$ such that:

```math
\rho = \rho(\boldsymbol{\theta}) = \frac{\mathbb{1}}{d} + \sum_{n=1}^{d^2-1} \theta_n \omega_n \qquad (1)
```

where $\omega_n$ denotes one of the generalized Gell-Mann matrices. The vector $\boldsymbol{\theta}$ is called the generalized Bloch vector. The Bloch vector can be calculated from the matrix $\rho$ using the function [`traceless_vectorization`](@ref). Reciprocally, the density matrix can be reconstructed from the Bloch vector using [`density_matrix_reconstruction`](@ref).

## Linear Inversion estimator

By applying $\mathrm{Tr} (\Pi_m \bullet)$ to both sides of equation (1), we arrive at the linear system:

```math
\mathbf{q} = T \boldsymbol{\theta}
```

where $q_m = p_m - \mathrm{Tr} \Pi_m / d$ and $T$ is a matrix with entries $T_{mn} = \mathrm{Tr} \left(\Pi_m \omega_n \right)$.

If the matrix $T$ is injective, the linear system has a unique solution given by:

```math
\boldsymbol{\theta} = (T^\dagger T)^{-1} T^\dagger\boldsymbol{q} \qquad (2)
```

In this case, the POVM is said to be **informationally complete**. Otherwise, if $T$ is not injective, the solution is no longer unique and the POVM is said to be informationally incomplete.

This defines the simplest measurement type in our package: given an iterable representing a measurement (each element of the iterable must be a matrix representing the POVM element, or a vector. In the latter case, the vector is interpreted as a projector). Then one calls [`assemble_measurement_matrix`](@ref) to construct the matrix $T$. This is passed to the functions that perform the tomography.

This would then be the simplest application of our package:
```julia
itr = ... # an iterable of matrices or vectors representing the POVM elements
probs = ... # a vector of probabilities corresponding to the outcomes of the measurement

μ = assemble_measurement_matrix(itr) # assemble the measurement matrix T
ρ, θ = estimate_state(probs, μ, LinearInversion()) # estimate the state (both the density matrix and the Bloch vector)
```

It is important to observe that the probabilities in the inversion formula (2) are not directly measurable and must be estimated. The simplest method substitutes them with observed experimental frequencies $\hat{p}_m = N_m / N$, where $N_m$ is the number of times outcome $m$ was observed and $N = \sum_m N_m$ is the total number of observations. This estimation of the probabilities is inherently noisy, due to finite number of observations or experimental imperfections, which can lead to inaccuracies in the reconstructed state. Most importantly, it might be the case that the state reconstructed from the estimated probabilities is not a valid density matrix, i.e., it may not be positive semi-definite. Therefore, the [`estimate_state`](@ref) function will apply an algorithm to project the reconstructed state onto the closest valid density matrix [smolin_efficient_2012](@cite), ensuring that the output is always a valid quantum state.

## Other estimation Methods

### Variants of Linear Inversion

While the basic linear inversion approach provides a straightforward solution to quantum state tomography, different computational strategies can be employed depending on the specific requirements of the problem. The package implements three variants of linear inversion, each with distinct advantages and trade-offs.

#### LinearInversion

The standard [`LinearInversion`](@ref) method solves the linear system $T \boldsymbol{\theta} = \boldsymbol{q}$ directly using QR decomposition via Julia's [`\`](https://docs.julialang.org/en/v1/base/math/) operator. This is the most straightforward implementation that solves the system without explicitly computing matrix inverses.

**Theory**: The method solves the linear system $T \boldsymbol{\theta} = \boldsymbol{q}$ directly using QR decomposition, which provides a numerically stable least-squares solution for overdetermined systems.

**Motivation**: This approach is ideal for single-shot tomography or when the measurement setup changes between different state estimations. It requires no preprocessing and is numerically stable due to the QR decomposition.

**Pros**: 
- Simple and robust implementation
- No memory overhead for preprocessing
- Numerically stable via QR decomposition
- Suitable for one-time or infrequent state estimations

**Cons**:
- Computational overhead for repeated estimations with the same measurement
- Slower than optimized variants for multiple estimations

#### PreAllocatedLinearInversion

The [`PreAllocatedLinearInversion`](@ref) method precomputes the pseudoinverse $T^+$ of the measurement matrix during initialization, storing it for subsequent state estimations. This eliminates the need to solve the linear system repeatedly.

**Theory**: By precomputing $T^+$, the state estimation reduces to a simple matrix-vector multiplication: $\boldsymbol{\theta} = T^+ \boldsymbol{q}$, where $\boldsymbol{q}$ is adjusted for the trace part.

**Motivation**: This variant is designed for scenarios where many state estimations are performed using the same measurement setup, such as real-time tomography.

**Pros**:
- Significant speedup for repeated estimations
- Predictable computational cost per estimation
- Optimal for batch processing or real-time applications

**Cons**:
- Memory overhead for storing the precomputed pseudoinverse
- Requires reinitialization if the measurement setup changes
- Initial setup cost for computing the pseudoinverse

#### NormalEquations

The [`NormalEquations`](@ref) method solves the normal equations $T^\dagger T \boldsymbol{\theta} = T^\dagger \boldsymbol{q}$ directly, which can be computationally advantageous when the number of measurements is much larger than the number of parameters.

**Theory**: Instead of computing the pseudoinverse, this method forms the normal equations explicitly: $(T^\dagger T) \boldsymbol{\theta} = T^\dagger \boldsymbol{q}$. The system matrix is square and typically smaller than the original overdetermined system.

**Motivation**: For measurements with a large number of POVM elements (large $M$) relative to the Hilbert space dimension, solving the normal equations can be more efficient than direct pseudoinverse computation.

**Pros**:
- Efficient for overcomplete measurements (large $M$ relative to $d^2-1$)
- Reduced computational complexity in high-redundancy scenarios
- Lower memory requirements during computation

**Cons**:
- Potentially less numerically stable than QR-based methods
- Can amplify numerical errors if $T^\dagger T$ is ill-conditioned
- May not provide advantages for moderately sized problems

#### Choosing the Right Variant

The choice between these variants depends on the specific use case:

- Use **LinearInversion** for single estimations or when numerical stability is paramount
- Use **PreAllocatedLinearInversion** when performing many estimations with the same measurement setup
- Use **NormalEquations** when dealing with highly redundant measurements where $M \gg d^2-1$

### Maximum Likelihood Estimator

For situations with limited observations (e.g., photocount regime), the maximum likelihood estimator is more appropriate:

```math
\hat{\rho}_{ML} = \arg\max \mathcal{L}(N_1, \ldots, N_M | \rho)
```

where the likelihood function is:

```math
\mathcal{L}(N_1, \ldots, N_M | \rho) = \prod_{m=1}^M \left[\mathrm{Tr} \left(\Pi_m \rho \right)\right]^{N_m}
```

This estimator ensures that the output is always a valid density operator and is more robust for small sample sizes.

### Bayesian Inference