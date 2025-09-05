#
# Monta a matriz global de rigidez
#
function RigidezGlobal(nn,ne,MAT,ESP,XY,IJ,ρ,P)

    # Inicializa a matriz de rigidez
    K = spzeros(2*nn,2*nn)

    # Loop pelos elementos finitos
    for e = 1:ne

        #  Recupera as propriedades do material do elemento 
        E = MAT[e,1]
        ν = MAT[e,2]

        # Espessura 
        te = ESP[e]
        
        # Coordenadas dos nós do elemento 
        X,Y = MontaXY(e,IJ,XY) 

        # Matriz de rigidez local do elemento
        Ke = MontaKe(X,Y,E,ν,te)

        # Determina os gls globais do elemento
        ge = Montagls(e,IJ)

        # Contribuição da variável de projeto
        mult = ρ[e]^P

        # Sobreposição da matriz de rigidez
        K[ge,ge] .+= Ke*mult

    end #e

    # Retorna a matriz de rigidez global 
    return K

end

#
# Vetor de forças global
#
function ForcaCglobal(nn,nf,FC)

    # Aloca o vetor de forças concentradas
    F = zeros(2*nn)

    # Laço pelas linhas de FC
    for i=1:nf

        # Nó
        no = FC[i,1]

        # gl local (1/X e 2/Y)
        gll = FC[i,2]

        # gl global
        glg = Int(2(no-1)+gll)

        # Valor da força
        valor = FC[i,3]

        # Sobrepõe no vetor global
        F[glg] = F[glg] + valor

    end

    # Retorna o vetor de foças concentradas 
    return F

end