module QuantumMeasurements

using LinearAlgebra, Random, OnlineStats, Distributions
import Base.rand
import LinearAlgebra.isposdef!
using Base.Threads: @spawn

include("gell_mann_matrices.jl")
export GellMannMatrices

include("vectorization.jl")
export traceless_vectorization!, traceless_vectorization, traceless_reconstruction!, traceless_reconstruction,
    density_matrix_reconstruction!, density_matrix_reconstruction,
    vectorization!, vectorization, reconstruction!, reconstruction

include("measurement.jl")
export get_traceless_part, get_trace_part, get_dim, get_num_outcomes, get_measurement_type,
    assemble_measurement_matrix, get_probabilities!, get_probabilities, fisher, empty_measurement,
    update_measurement!, multithreaded_update_measurement!, filter_measurement

include("proportional_measurement.jl")
export ProportionalMeasurement

include("linear_inversion.jl")
include("max_likelihood.jl")
include("bayesian_inference.jl")
export estimate_state, LinearInversion, PreAllocatedLinearInversion, NormalEquations,
    MaximumLikelihood, BayesianInference

include("misc.jl")
export fidelity, project2density!, project2density, project2pure, polarization_state,
    get_noisy_probabilities, simulate_outcomes

include("random_matrices.jl")
export HaarUnitary, HaarVector, Simplex, ProductMeasure, GinibreEnsamble

using PrecompileTools: @setup_workload, @compile_workload
include("precompile.jl")

end
