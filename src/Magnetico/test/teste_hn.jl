@testset "Quad4" begin

    # Chama o programa para o problema de teste
    ϕ,H = Analise("malhas/1D_hn_quad.msh",output=false)

    # Para esse problema, o valor mínimo de ϕ deve ser 
    @test minimum(ϕ) ≈ -2.0

    # e o valor de Hx deve ser 
    @test minimum(H[:,1]) ≈ 1.0

    # e o valor de Hy deve ser 
    @test isapprox(maximum(abs.(H[:,2])), 0.0, atol=1E-6)


end

@testset "Tri3" begin

    # Chama o programa para o problema de teste
    ϕ,H = Analise("malhas/1D_hn_tri.msh",output=false)

    # Para esse problema, o valor mínimo de ϕ deve ser 
    @test minimum(ϕ) ≈ -2.0

    # e o valor de Hx deve ser 
    @test minimum(H[:,1]) ≈ 1.0

    # e o valor de Hy deve ser 
    @test isapprox(maximum(abs.(H[:,2])), 0.0, atol=1E-6)


end