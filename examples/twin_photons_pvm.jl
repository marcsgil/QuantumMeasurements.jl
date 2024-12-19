# # Twin Photons
# In this example, we perform the tomography of the polarization state of a pair of twin photons.
# The photons are generated from SPDC and are entangled in polarization.

# We will first load up a few packages. 
# We use [CSV.jl](https://github.com/JuliaData/CSV.jl) (you might need to install it) 
# in order to read the real experimental data.

using CSV, LinearAlgebra, QuantumMeasurements

# We then read the [experimental data](@__REPO_ROOT_URL__/examples/counts.csv), 
# which is located in the [`examples`](@__REPO_ROOT_URL__/examples) directory of the repository.

dir = pkgdir(QuantumMeasurements);
path = joinpath(dir, "examples", "counts.csv");
file = CSV.File(path, header=false);

# The following code translates the contents of the file to a form that will be useful to us.

# This is just a function to parse a string into a complex number
parse_c(s) = parse(ComplexF32, s);

# The outcomes are contained in row 4, and we normalize to get frequencies
outcomes = [Float32(parse_c(row[4])) for row in file]
freqs = normalize(outcomes, 1);

# These are the coefficients defining the measurement, which are coincidence counts after the beam gets 
# projected into a pair of polarization states

# These vectors define the polarization projector for each of the two detectors
ψ1 = [[parse_c(row[5]), parse_c(row[6])] for row ∈ file]
ψ2 = [[parse_c(row[7]), parse_c(row[8])] for row ∈ file];

# We now assemble the projectors on the two-photon state using [`kron`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#Base.kron).
# We divide by 3 to ensure that the sum of the PVM projectors is the identity matrix.
pvm = [kron(pair...) / 3 for pair in zip(ψ1, ψ2)];

# The measurement matrix is assembled:
μ = assemble_measurement_matrix(pvm);

# Finally, we choose the [`MaximumLikelihood`](@ref) method and make our estimation:
method = MaximumLikelihood()
ρ_pred = estimate_state(freqs, μ, method)[1]

# We can see that we are very close to the Bell state, which was the desired outcome of the experiment:
ψ = [1 + 0im, 0, 0, 1] / √2
fidelity(ψ, ρ_pred) ≥ 0.99