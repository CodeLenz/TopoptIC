@testset "Quad4" begin

    # Chama o programa para o problema de teste
    ϕ,H = Analise("malhas/1D_rhom_quad.msh",output=false)

    # Para esse problema, o valor mínimo de ϕ deve ser 
    @test maximum(ϕ) ≈ 2.0

    # e o valor de Hx deve ser 
    @test isapprox(minimum(H[:,1]),-2.0,atol=1E-2)

    # e o valor de Hy deve ser 
    @test isapprox(maximum(abs.(H[:,2])), 0.0, atol=1E-2)

end

@testset "Tri3" begin

    # Chama o programa para o problema de teste
    ϕ,H = Analise("malhas/1D_rhom_tri.msh",output=false)

    # Para esse problema, o valor mínimo de ϕ deve ser 
    @test isapprox(maximum(ϕ), 2.0, atol=1E-3)

    # e o valor de Hx deve ser 
    @test isapprox(minimum(H[:,1]),-2.0,atol=1E-2)

    # e o valor de Hy deve ser 
    @test isapprox(maximum(abs.(H[:,2])), 0.0, atol=1E-2)

end