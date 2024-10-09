pol_states = (polarization_state(s) for s ∈ [:H, :V, :D, :A, :R, :L])

μ_pvm = assemble_measurement_matrix(pol_states)
μ_povm = assemble_measurement_matrix(map(x -> x * x', pol_states))

@testset "POVM x PVM" begin
    @test μ_pvm ≈ μ_povm
end