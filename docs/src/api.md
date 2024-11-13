# API

## Measurements

### Measurement interface

The methods in this section forms the interface for measurements. The user can implement their own measurements by following the interface.

```@autodocs
Modules = [QuantumMeasurements]
Pages = ["measurement.jl"]
Filter = f -> f != assemble_measurement_matrix
```

If [`get_measurement_type`](@ref), [`get_num_outcomes`](@ref) and [`get_probabilities!`](@ref) are implemented, one has also access to [`get_probabilities`](@ref):

```@docs
get_probabilities
```

### Measurement matrix

The simplest measurement type is just a matrix:

```@docs
assemble_measurement_matrix(itr)
```

### Proportional Measurements

```@docs
ProportionalMeasurement
```

## Tomography Methods

This package contains several methods for quantum state tomography. They should be passed to the following function in order to perform the tomography:

```@docs
estimate_state
```

The available methods are:

```@docs
LinearInversion
PreAllocatedLinearInversion
NormalEquations
MaximumLikelihood
BayesianInference
```

## Gell-Mann Matrices

```@docs
GellMannMatrices
```

## Random

```@autodocs
Modules = [QuantumMeasurements]
Pages = ["random_matrices.jl"]
```

## Misc

```@autodocs
Modules = [QuantumMeasurements]
Pages = ["misc.jl"]
```