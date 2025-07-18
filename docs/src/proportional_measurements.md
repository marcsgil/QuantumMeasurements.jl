# [Proportional Measurements](@id prop_meas)

## Obstructed Measurements

### Conditional Probabilities

When measurements are obstructed, the experimental setup only allows estimation of conditional probabilities:

```math
p_{m | D} = \frac{p_m}{p_D} = \frac{\mathrm{Tr} \left[\rho \Pi(\mathcal{R}_m) \right]}{\mathrm{Tr} \left(\rho g \right)}
```

where $g = \sum_{m=1}^{M} \Pi(\mathcal{R}_m)$ is the total detection operator and $p_D = \mathrm{Tr}(\rho g)$ is the total detection probability.

### State Transformation

The measurement process transforms the state through:

```math
\rho \mapsto \tilde{\rho} = \frac{A \rho A^\dagger}{\mathrm{Tr} \left(\rho g\right)}
```

where $A$ satisfies $A^\dagger A = g$. The original state can be recovered by:

```math
\rho = \frac{A^{-1} \tilde{\rho} (A^{-1})^\dagger}{\mathrm{Tr} \left[A^{-1} \tilde{\rho} (A^{-1})^\dagger \right]}
```

This framework allows tomographic reconstruction even in the presence of significant obstructions, extending the applicability of the method to practical scenarios such as free-space optical communication.