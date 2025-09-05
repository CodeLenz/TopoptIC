function Derivada_C(ne,ρ,P,U,MAT,ESP,XY,IJ)

    # Inicializa do vetor com as derivadas
    dC = zeros(ne)

    # Loop pelos elementos finitos
    for m = 1:ne

        #  Recupera as propriedades do material do elemento 
        E = MAT[m,1]
        ν = MAT[m,2]

        # Espessura 
        te = ESP[m]
        
        # Coordenadas dos nós do elemento 
        X,Y = MontaXY(m,IJ,XY) 

        # Matriz de rigidez local do elemento
        Ke0 = MontaKe(X,Y,E,ν,te)

        # Determina os gls globais do elemento
        ge = Montagls(m,IJ)

        # Contribuição da variável de projeto
        mult = P*ρ[m]^(P-1)

        # Extrai o deslocamento do elemento m
        um = U[ge]

        # Calcula a derivada
        dC[m] = -transpose(um)*(mult*Ke0)*um

    end #e

    # Retorna a derivada
    return dC

end