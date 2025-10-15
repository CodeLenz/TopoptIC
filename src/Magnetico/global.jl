#
# Monta a matriz global K
#
function Monta_K(ne,coord,connect,materials)  

    # Aloca vetores para a montagem eficiente 
    # das matrizes esparsas
    I = Int64[]
    J = Int64[]
    VK = Float64[]
    
    # Loop pelos elementos
    for ele=1:ne

        # Material
        mat = connect[ele,2] 

        # Find material μ
        μ = materials[mat,1]
        
        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X, Y e Z para este elemento
        nos, X = Nos_Coordenadas(ele,et,coord,connect) 

        # Monta as matrizes dos elementos
        Ke = Ke_bi4(μ,X)

        #=
        if et==3
           Ke, Me = KMe_bi4(iρ,iκ,X)
        elseif et==2
           Ke, Me = KMe_tri3(iρ,iκ,X)
        elseif et==4
            Ke, Me = KMe_tet4(iρ,iκ,X)   
        elseif et==5
           Ke, Me = KMe_hex8(iρ,iκ,X)
        elseif et==7
            Ke, Me = KMe_pyr5(iρ,iκ,X)    
        else
            error("Elemento não definido")
        end
        =#
 
        # Sobreposição das locais nas globais
        for i in LinearIndices(nos)
            ni = nos[i]
            for j in LinearIndices(nos)
                push!(I,ni)
                push!(J,nos[j])
                push!(VK,Ke[i,j])
            end
        end

    end

    # Retorna as matrizes globais
    return sparse(I,J,VK)

end
