@testset "Matriz_N" begin

    # Testa as funções N nas posições dos nós 
    #
    # N(-1,-1) = [1 0 0 0] 
    #
    #
    N =  Matriz_N_bi4(-1,-1)
    @test all(N .≈ [1 0 0 0])

    # Testa as funções N nas posições dos nós 
    #
    # N(1,-1) = [0 1 0 0] 
    #
    #
    N =  Matriz_N_bi4(1,-1)
    @test all(N .≈ [0 1 0 0])

    # Testa as funções N nas posições dos nós 
    #
    # N(1,1) = [0 0 1 0] 
    #
    #
    N =  Matriz_N_bi4(1,1)
    @test all(N .≈ [0 0 1 0])

    # Testa as funções N nas posições dos nós 
    #
    # N(-1,1) = [0 0 0 1] 
    #
    #
    N =  Matriz_N_bi4(-1,1)
    @test all(N .≈ [0 0 0 1])

end
