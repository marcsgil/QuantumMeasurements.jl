# QuantumMeasurements.jl

[![DOI](https://zenodo.org/badge/870243634.svg)](https://doi.org/10.5281/zenodo.14009331)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcsgil.github.io/QuantumMeasurements.jl/)

**A comprehensive Julia package for quantum state tomography with multiple estimation methods.**

QuantumMeasurements.jl provides a complete toolkit for reconstructing quantum states from experimental measurement data. The package implements state-of-the-art algorithms ranging from fast linear inversion to sophisticated Bayesian inference, making it suitable for both routine tomography and cutting-edge research applications.

## ✨ Key Features

- **Multiple Estimation Methods**: Linear inversion, maximum likelihood, and Bayesian inference
- **Flexible Measurements**: Support for standard POVMs and proportional measurements
- **High Performance**: Optimized algorithms with pre-allocated variants for repeated measurements  
- **Uncertainty Quantification**: Full posterior distributions with error bars via MCMC
- **Experimental Validation**: Handles obstructed measurements and incomplete data
- **Testing Utilities**: Built-in random state generators for validation and benchmarking
- **Comprehensive Documentation**: Detailed theory, examples, and practical guides

## 🚀 Quick Start

### Installation

```julia
using Pkg
Pkg.add("QuantumMeasurements")
```

### Basic Usage

```julia
using QuantumMeasurements

# Define measurement basis (polarization tomography example)
basis = [polarization_state(s) / √3 for s ∈ (:H, :V, :D, :A, :R, :L)]
measurement = assemble_measurement_matrix(basis)

# Experimental outcomes (photon counts)
outcomes = [100, 0, 50, 50, 50, 50]

# Perform tomography
method = LinearInversion()
ρ_estimated = estimate_state(outcomes, measurement, method)[1]
```

## 📊 Estimation Methods

| Method | Best For | Speed | Uncertainty |
|--------|----------|-------|-------------|
| `LinearInversion` | High photon counts | Fast | No |
| `MaximumLikelihood` | Low photon counts | Medium | No |
| `BayesianInference` | Full uncertainty analysis | Slow | Yes |

## 📖 Documentation

- **[Documentation Home](https://marcsgil.github.io/QuantumMeasurements.jl/)**: Complete documentation with tutorials and examples
- **[Quick Start Guide](https://marcsgil.github.io/QuantumMeasurements.jl/quick_start/)**: Get started with basic tomography
- **[Method Selection](https://marcsgil.github.io/QuantumMeasurements.jl/choosing_methods/)**: Choose the right algorithm for your experiment
- **[Examples](https://marcsgil.github.io/QuantumMeasurements.jl/examples/)**: Practical applications and use cases

## 🧪 Example Applications

- **Polarization State Tomography**: Single and two-photon systems
- **Spatial Light Tomography**: Orbital angular momentum states with obstructions
- **Photon Counting Experiments**: Low-count regime with proper statistical treatment
- **Real-time Tomography**: High-throughput applications with optimized algorithms

## 🔬 Research Applications

This package has been used in research on:
- Spatial structure tomography of light with physical obstructions
- Informationally complete measurements using astigmatic transformations

## 📝 Citation

If you use QuantumMeasurements.jl in your research, please cite:

```bibtex
@software{marcsgil_2024_14009332,
  author       = {Gil de Oliveira, Marcos},
  title        = {QuantumMeasurements.jl},
  month        = oct,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {0.1.0},
  doi          = {10.5281/zenodo.14009332},
  url          = {https://github.com/marcsgil/QuantumMeasurements.jl}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
