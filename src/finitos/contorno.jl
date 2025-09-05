#
# Mapeamento dos nós locais das faces do elemento
#
function FACES(face)

    if face==1
       p1 = 1
       p2 = 2
    elseif face==2
       p1 = 2
       p2 = 3
    elseif face==3
       p1 = 3
       p2 = 4
    elseif face==4
       p1 = 4
       p2 = 1
    else
       error("FACES:: face $face não existe")
       #error("FACES:: face ",face," não existe")
    end   

    # retorna os indicadores dos nós das faces
    return p1,p2

end

#
# Monta a matriz N para integração nas faces do elemento 
# bilinear isop de 4 nós
#
function MontaN(ζ,p1,p2)

    # Aloca a matriz de saída
    N = zeros(2,8)

    # Funções de interpolação não nulas
    N1 = (1/2)*(1-ζ)
    N2 = (1/2)*(1+ζ)

    # Posiciona na matriz N
    N[1,2*p1-1] = N1
    N[2,2*p1]   = N1
    N[1,2*p2-1] = N2
    N[2,2*p2]   = N2
    
    # Return N
    return N

end

#
# Determina o sistema  local v × n de uma face 
# do elemento
#
function Facesnv(p1,p2,e,IJ,XY)

    # Recupera os nós (global) da face
    n1 = IJ[e,p1]
    n2 = IJ[e,p2]

    # Recupera as coordenadas desses nós
    x1 = XY[n1,1]
    y1 = XY[n1,2]
    x2 = XY[n2,1]
    y2 = XY[n2,2]
    
    # Calcula dos deltas
    ΔX = x2 - x1
    ΔY = y2 - y1

    # Comprimento da face (no sistema original)
    L = sqrt(ΔX^2 + ΔY^2)

    # Determinante do Jacobiano 
    dJ = L/2

    # Vetor tangente à face
    v = (1/L)*[ΔX ; ΔY]

    # Vetor normal à face
    n = (1/L)*[ΔY;-ΔX]

    # Retorna dJ, n, v
    return dJ, n, v

end

#
# Integral de contorno em uma das faces do elemento
#
# carregamento pode ser tanto p (normal à face) quanto um τ
# (tangencial à face). 
#
# vetor_dir_local pode ser n (p) ou v (τ)
#
function Ftblinear(dJ,vetor_dir_local,p1,p2,carregamento,te)

    # Pontos da quadraduta de Gauss-Legendre
    pg = (1/sqrt(3))*[-1;1]

    # Aloca o vetor de saída
    F = zeros(8)

    # Laço de integração da quadratura
    for i=1:2

        # Recupera o ponto de Gauss
        ζ = pg[i]

        # Monta a matriz N neste ponto/face
        N = MontaN(ζ,p1,p2)

        # Acumula a integral
        F .= F .+ carregamento*transpose(N)*vetor_dir_local*dJ*te

    end #i

    # Retorna o vetor de forças na face
    return F

end

#
# Monta o vetor global de forças devido as forças aplicadas 
# nas faces dos elementos
#
function ForcatGlobal(nn,ESP,XY,IJ,np,P)

    # Inicializa o vetor global 
    F = zeros(2*nn)

    # Laço pelas linhas de P
    for i=1:np

        # Recupera as informações da linha de P
        ele  = Int(P[i,1])
        face = Int(P[i,2])
        dir  = Int(P[i,3])
        val  = P[i,4]

        # Recupera a espessura do elemento 
        te = ESP[ele]

        # Ponteiros da face
        p1,p2 = FACES(face)

        # Vetores da face
        dJ,n,v = Facesnv(p1,p2,ele,IJ,XY)

        # Determina se passamos n ou v para  Fte
        if dir==1
           vetor = n
        else
           vetor = v
        end

        # Calcula a integral na face
        Fte = Ftblinear(dJ,vetor,p1,p2,val,te)

        # Determina os graus de liberdade globais do elemento 
        ge = Montagls(ele,IJ)

        # Posiciona o vetor local 8 × 1 no vetor global 
        for j=1:8
            F[ge[j]] = F[ge[j]] + Fte[j]
        end
        

    end #i

    # Retorna o vetor global de forças distribuídas
    return F

end