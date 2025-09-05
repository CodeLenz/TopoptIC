#
# Derivada das funções de interpolação 
# em relação a r e s
#
function MontadN(r,s)

    # Matriz com as derivadas das funções de interpo
    # primeira linha em relação a r
    # segunda linha em relação a s
    dNrs = (1/4)*[s-1   1-s   1+s -(1+s) ;
                  r-1 -(1+r)  1+r   1-r ]

    # Retorna a matriz
    return SMatrix{2,4}(dNrs)

end


#
# Monta a matriz J em um ponto (r,s)
#
function MontaJ(dNrs,X,Y)

    # Aloca a matriz J
    J = zeros(2,2)

    # Loop do somatório
    for i=1:4

        # dx/dr
        J[1,1] = J[1,1] + dNrs[1,i]*X[i]

        # dy/dr
        J[1,2] = J[1,2] + dNrs[1,i]*Y[i]

        # dx/ds
        J[2,1] = J[2,1] + dNrs[2,i]*X[i]

        # dy/ds
        J[2,2] = J[2,2] + dNrs[2,i]*Y[i]

    end #i

    # Retorna a matriz J no ponto (r,s)
    return SMatrix{2,2}(J)

end

#
# Corrige as derivadas dNrs para dNxy em um ponto 
# (r,s) do elemento
#
function CorrigedN(dNrs,J)

    # Inverte a matriz
    invJ = inv(J)

    # Calcula a correção 
    dNxy = invJ*dNrs

    # Retorna a correção 
    return dNxy

end

#
# Monta a matriz B no ponto (r,s)
#
function MontaB(dNxy)

    # Aloca a matriz B
    B = zeros(3,8)

    # Indicador de começo de bloco
    c = 1

    # Loop pelas funções de interpolação 
    for i=1:4
 
       # Derivada em relação a x
       B[1,c] = dNxy[1,i]

       # Derivada em relação a y
       B[2,c+1] = dNxy[2,i]

       # Derivada em relação a y
       B[3,c] = dNxy[2,i]

       # Derivada em relação a x
       B[3,c+1] = dNxy[1,i]

       # Atualiza o início do bloco
       c = c + 2

    end

    # Retorna a matriz B no ponto (r,s)
    return SMatrix{3,8}(B)

end



#
# Monta XY
#
function MontaXY(e,IJ,XY)

    # Aloca os vetores de saída
    X = zeros(4)
    Y = zeros(4)

    # Loop pelos nós deste elemento
    for n=1:4

        # Descobre quem é o nó i do elemento 
        no = IJ[e,n]

        # Recupera as coordeandas deste nó 
        # na matriz XY
        X[n] = XY[no,1]
        Y[n] = XY[no,2]

    end

    # Retorna os vetores 
    return X,Y

end


#
# Monta o vetor com os gls globais do elemento
#
function Montagls(e,IJ)

    # Aloca o vetor de saída
    gls = zeros(Int64,8)

    # Contador 
    c = 1

    # Loop pelos nós
    for i=1:4

        # Recupera o nó 
        no = IJ[e,i]

        # Loop pelos gls
        for j=1:2

            # Calcula o gl 
            gls[c] = 2*(no-1)+j

            # Incrementa o contador 
            c = c + 1

        end #j

    end #i

    # Retorna o vetor com os gls globais do elemento 
    return gls

end

#
# Monta a matriz de rigidez local do elemento
#
function MontaKe(X,Y,E,ν,te)

    # Aloca a matriz
    Ke = @MMatrix zeros(8,8)

    # Monta a matriz constitutiva
    Cv = MontaCEPT(E,ν)

    # Define os pontos de Gauss-Legendre
    pg = (1/sqrt(3))*[-1  1 1 -1 ;
                      -1 -1 1  1  ]

    # Laço pelos pontos de Gauss
    for i=1:4

       # Recupera as coordenadas do PG
       r = pg[1,i]
       s = pg[2,i]

       # Derivadas das funções de interpolação 
       # em relação a rs neste ponto
       dNrs = MontadN(r,s)

       # Calcula a matriz J no ponto
       J = MontaJ(dNrs,X,Y)

       # Calcula o determinante no ponto 
       dJ = det(J)

       # Teste de consistência de conectividades
       dJ > 0 || error("dJ negativo $dJ")

       # Mapeia as derivadas para xy
       dNxy = CorrigedN(dNrs,J)
  
       # Monta a matriz B no ponto
       B = MontaB(dNxy)

       # Acumula a matriz (somatório da integração)
       Ke .+=  transpose(B)*Cv*B*dJ*te

    end #i

    # Retorna a matriz de rigidez do elemento
    return Symmetric(Ke)
    
end


#
# Calcula a área de do elemento
#
function Area_elemento(e,IJ,XY)

    # Inicializa a área do elemento
    Ae = 0.0

    # Coordenadas dos nós do elemento
    X,Y = MontaXY(e,IJ,XY)

    # Define os pontos de Gauss-Legendre
    pg = (1/sqrt(3))*[-1  1 1 -1 ;
                      -1 -1 1  1  ]

    # Laço pelos pontos de Gauss
    for i=1:4

       # Recupera as coordenadas do PG
       r = pg[1,i]
       s = pg[2,i]

       # Derivadas das funções de interpolação 
       # em relação a rs neste ponto
       dNrs = MontadN(r,s)

       # Calcula a matriz J no ponto
       J = MontaJ(dNrs,X,Y)

       # Calcula o determinante no ponto 
       dJ = det(J)

       # Acumula a área
       Ae +=  dJ

    end #i

    # Retorna a área do elemento
    return Ae
    
end