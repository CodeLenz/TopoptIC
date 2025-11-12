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
        #Ke = Ke_bi4(μ,X)

        if et==3
           Ke = Ke_bi4(μ,X)
        elseif et==2
           Ke = Ke_tri3(μ,X)
        elseif et==4
            Ke = Ke_tet4(μ,X)
        elseif et==5
           Ke = Ke_hex8(μ,X)
        elseif et==7
            Ke = Ke_pyr5(μ,X)
        else
            error("Global::Elemento não definido")
        end
    
 
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
