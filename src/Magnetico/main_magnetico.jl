function Analise(meshfile::String)
    
    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, bn, _ = Parsemsh(meshfile)

    # Precisamos de um material
    if isempty(materials)
        error("Analise:: at least one material is necessary")
    end

    # Calcula a matriz de rigidez global
    K = Monta_K(ne,coord,connect,materials)

    # Posições que não precisam ser calculadas no sistema de equações 
    livres = setdiff(collect(1:nn),nodes_open)
    
    # Vetor de forças 
    P = Vetor_P(nn,bn,coord,connect)

    # Aloca um vetor Φ com a dimensão completa
    Φ = zeros(nn)

    # Soluciona o sistema de equações para os gls livres
    Φ[livres] .= -K[livres,livres]\P[livres]

   
end