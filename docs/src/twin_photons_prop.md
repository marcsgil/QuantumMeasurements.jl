```@meta
EditURL = "../../examples/twin_photons_prop.jl"
```

# Twin photons with proportional measurement

In this example, we perform a similar tomography as in [Twin Photons](@ref).
The difference is that, now, the measurement performed is not a PVM.
The data from this example was taken from [PhysRevA.64.052312](@cite)

We load the necessary packages and define our polarization states

````@example twin_photons_prop
using LinearAlgebra, QuantumMeasurements

h = polarization_state(:H)
v = polarization_state(:V)
d = polarization_state(:D)
a = polarization_state(:A)
r = polarization_state(:R)
l = polarization_state(:L);
nothing #hide
````

The main difference arises here, in the specification of our measurement:

````@example twin_photons_prop
projectors = [kron(h, h),
    kron(h, v),
    kron(v, v),
    kron(v, h),
    kron(r, h),
    kron(r, v),
    kron(d, v),
    kron(d, h),
    kron(d, r),
    kron(d, d),
    kron(r, d),
    kron(h, d),
    kron(v, d),
    kron(v, l),
    kron(h, l),
    kron(r, l),];
nothing #hide
````

Notice that the sum of the corresponding projectors is not even proportional to the identity matrix:

````@example twin_photons_prop
sum(x -> x * x', projectors)
````

Therefore, this measurement is not a PVM/POVM.

As discussed in [PhysRevA.64.052312](@cite), this occurred because of experimental ease:
This measurement only contains 16 settings (the minimum amount), as opposed to the
36 settings present in the [Twin Photons](@ref) example.
Notice, also, that only one of the polarization states changes between each successive measurement,
so that only one set of wave plates needs to be adjusted.

To deal with this kind of situation, we can build a [`ProportionalMeasurement`](@ref):

````@example twin_photons_prop
μ = ProportionalMeasurement(projectors);
nothing #hide
````

The theory behind this is better explained at [Proportional Measurements](@ref prop_meas).

From now on, everything follows the standard procedure.

We specify the experimental outcomes

````@example twin_photons_prop
outcomes = [34749, 324, 35805, 444, 16324, 17521, 13441,
    16901, 17932, 32028, 15132, 17238, 13171, 17170, 16722, 33586];
nothing #hide
````

We define our tomography method and make our prediction

````@example twin_photons_prop
method = MaximumLikelihood()
ρ_pred = estimate_state(outcomes, μ, method)[1]
````

Once again, we aim to be close to a Bell state:

````@example twin_photons_prop
ψ = [1 + 0im, 0, 0, 1] / √2
0.96 ≤ fidelity(ρ_pred, ψ) ≤ 0.97
````

```@bibliography
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

