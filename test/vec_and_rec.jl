dims = 10:10:100

@testset "Vectorization and reconstruction tests" begin
    for dim ∈ dims
        θs = Vector{Float32}(undef, dim^2 - 1)
        xs = Vector{Float32}(undef, dim^2)
        ωs = GellMannMatrices(dim)
        σ = Matrix{ComplexF32}(undef, dim, dim)
        X = rand(ComplexF32, dim, dim)
        hermitianpart!(X)

        traceless_vectorization!(θs, X)
        @test [real(tr(X * ω)) for ω in ωs] ≈ θs
        @test sum(prod, zip(θs, ωs)) .+ tr(X) * I(dim) ./ dim ≈ X

        traceless_reconstruction!(σ, θs)
        @test σ ≈ X - tr(X) * I(dim) ./ dim

        ρ = X / tr(X)
        traceless_vectorization!(θs, ρ)
        density_matrix_reconstruction!(σ, θs)
        @test σ ≈ ρ

        @test reconstruction(vectorization(X)) ≈ X
    end
end