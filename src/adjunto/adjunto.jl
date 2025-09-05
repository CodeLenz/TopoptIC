#
# Calcula a parcela da derivada que é dependente do 
# vetor adjunto λ
#
function Derivada_termo_λ(ne,ρ,expP,U,λ,MAT,ESP,XY,IJ)

    # Inicializa do vetor com as derivadas
    dλ = zeros(ne)

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
        Km0 = MontaKe(X,Y,E,ν,te)

        # Determina os gls globais do elemento
        ge = Montagls(m,IJ)

        # Contribuição da variável de projeto
        mult = expP*ρ[m]^(expP-1)

        # Extrai o deslocamento do elemento m
        um = U[ge]
        λm = λ[ge]

        # Calcula a derivada
        dλ[m] = transpose(λm)*(mult*Km0)*um

    end #e

    # Retorna a derivada
    return dλ

end