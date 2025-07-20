# QuantumMeasurements.jl Documentation

Welcome to the documentation for **QuantumMeasurements.jl**, a comprehensive Julia package for quantum state tomography. This documentation will guide you through everything from basic usage to advanced theoretical foundations.

## What is QuantumMeasurements.jl?

QuantumMeasurements.jl reconstructs quantum states from experimental measurement data using statistically rigorous methods.

**Key Capabilities:**
- **Multiple estimation algorithms** from fast linear methods to full Bayesian inference
- **Flexible measurement handling** including incomplete and obstructed measurements
- **Uncertainty quantification** with proper error bars and confidence intervals

## Quick Navigation

### 🚀 **Getting Started**
New to quantum state tomography or the package? Start here:
- **[Installation](#installation)**: Get the package running
- **[Quick Start Guide](quick_start.md)**: Basic polarization tomography example
- **[Examples](examples/twin_photons_pvm.md)**: Practical applications

### 📚 **Theory and Methods**  
Understand the mathematics and choose the right approach:
- **[Mathematical Foundations](mathematical_foundations.md)**: Core theory and linear methods
- **[Maximum Likelihood Estimation](maximum_likelihood.md)**: Optimal point estimation
- **[Bayesian Inference](bayesian_inference.md)**: Full uncertainty quantification
- **[Proportional Measurements](proportional_measurements.md)**: Handling incomplete measurements

### 🛠 **Practical Guides**
Application-focused resources for real experiments:
- **[Choosing Estimation Methods](choosing_methods.md)**: Decision guide for your experimental setup
- **[Random State Generation](random_states.md)**: Testing and validation utilities

### 📖 **Reference**
Complete technical reference:
- **[API Documentation](api.md)**: All functions, types, and parameters
- **[References](references.md)**: Scientific citations and further reading

## Installation

QuantumMeasurements.jl is registered in the Julia General Registry. Install it using:

```julia
using Pkg
Pkg.add("QuantumMeasurements")
```

## Supported Julia Versions

This package supports Julia 1.10 and later versions.

## First Steps

Once installed, try the [Quick Start Guide](quick_start.md) for a hands-on introduction, or explore the [Examples](examples/twin_photons_pvm.md) to see the package in action with real experimental scenarios.

For questions about which method to use for your specific application, the [Choosing Estimation Methods](choosing_methods.md) guide provides practical decision-making guidance.