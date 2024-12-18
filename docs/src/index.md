```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: QuantumMeasurements
  tagline: A Julia package for quantum tomography
  actions:
    - theme: brand
      text: Quick Start 🚀
      link: quick_start
    - theme: alt
      text: API Reference 📚
      link: api
    - theme: alt
      text: View on GitHub
      link: https://github.com/LuxDL/Lux.jl
  image:
    src: /assets/logo.svg
    alt: QuantumMeasurements.jl

features:
  - icon: 🚀
    title: Fast & Extendible
    details: Lux.jl is written in Julia itself, making it extremely extendible. CUDA and AMDGPU are supported first-class, with experimental support for Metal and Intel GPUs.
    link: /introduction

  - icon: 🐎
    title: Powered by the XLA Compiler
    details: Lux.jl seamlessly integrates with Reactant.jl, to compiler models to run on CPU, GPU, TPU, and more.
    link: /manual/compiling_lux_models

  - icon: 🧑‍🔬
    title: SciML ❤️ Lux
    details: Lux is the default choice for all SciML packages, including DiffEqFlux.jl, NeuralPDE.jl, and more.
    link: https://sciml.ai/

  - icon: 🧩
    title: Uniquely Composable
    details: Lux.jl natively supports Arbitrary Parameter Types, making it uniquely composable with other Julia packages (and even Non-Julia packages).
    link: /api/Lux/contrib#Training
---
```

## How to Install QuantumMeasurements.jl?

Its easy to install Lux.jl. Since Lux.jl is registered in the Julia General registry,
you can simply run the following command in the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("Lux")
```

If you want to use the latest unreleased version of Lux.jl, you can run the following
command: (in most cases the released version will be same as the version on github)

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/LuxDL/Lux.jl")
```