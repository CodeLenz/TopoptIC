
#
# Monta a matriz N em um ponto (r,s)
#
function MontaN(r,s)

   # Calcula as funções no ponto 
   N1 = (1/4)*(1-r)*(1-s)
   N2 = (1/4)*(1+r)*(1-s)
   N3 = (1/4)*(1+r)*(1+s)
   N4 = (1/4)*(1-r)*(1+s)

   # Devolve a matriz
   [N1  0   N2  0  N3  0 N4 0 ;
     0  N1  0   N2 0 N3 0  N4]

end

#
# Monta o vetor local de forças de corpo devido a um vetor 
# b 2 × 1 [N/m³]
#
function MontaFBe(X,Y,te,b::Vector)

    # Aloca o vetor
    Fb = zeros(8)

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
  
       # Monta a matriz N no ponto
       N = MontaN(r,s)

       # Acumula o vetor 
       Fb .= Fb + transpose(N)*b*dJ*te

    end #i

    # Retorna o vetor
    return Fb
    
end


#
# Monta o vetor global de forças devido as forças de corpo
#
function ForcabGlobal(nn,ESP,XY,IJ,nfb,FB)

    # Inicializa o vetor global 
    F = zeros(2*nn)

    # Vetor b local com dimensão 2×1
    b = zeros(2) 

    # Laço pelas linhas de FB
    for i=1:nfb

        # Recupera as informações da linha de FB
        ele  = Int(FB[i,1])
        dir  = Int(FB[i,2])
        val  = FB[i,3]
 
        # vetor b local 2×1
        fill!(b,0)
        b[dir] = val
        
        # Coordenadas dos nós do elemento 
        X,Y = MontaXY(ele,IJ,XY) 

        # Recupera a espessura do elemento 
        te = ESP[ele]

        # Calcula a integral no domínio do elemento, retornando 
        # um vetor 8 × 1
        Fb = MontaFBe(X,Y,te,b)

        # Determina os graus de liberdade globais do elemento 
        ge = Montagls(ele,IJ)

        # Posiciona o vetor local 8 × 1 no vetor global 
        for j=1:8
            F[ge[j]] = F[ge[j]] + Fb[j]
        end
        
    end #i

    # Retorna o vetor global de forças de corpo
    return F

end