#
# Parametrização da tensão
#
function fσ(ρe, expP, expQ)

    ρe^(expP-expQ)

end

#
# Derivada da parametrização da tensão em relação a ρ
#
function dfσ(ρe, expP, expQ)

    (expP-expQ)*ρe^(expP-expQ-1)

end


#
# Calcula a restrição de tensão em um elemento 
#
function gσ(σeq,σY)

   return σeq/σY - 1

end

#
# Calcula a tensão equivalente em um elemento 
#
function Tensao_eq(σ, ϵ=1E-6)

    # Define a matriz de Voigt
    V = SMatrix{3,3}( [1.0 -0.5 0.0;
                      -0.5  1.0 0.0;
                       0.0  0.0 3.0] )

    # Calcula a tensão equivalente, evitando divisão por zero
    sqrt(transpose(σ)*V*σ + ϵ^2)

end


#
# Calcula a tensão no centro de um elemento finito 
#
function Tensao0_centro_elemento(ele,MAT,IJ,XY,U)

    #  Recupera as propriedades do material do elemento 
    E = MAT[ele,1]
    ν = MAT[ele,2]

    # Monta a matriz constitutiva
    Cv = MontaCEPT(E,ν)

    # Coordenadas dos nós do elemento 
    X,Y = MontaXY(ele,IJ,XY) 

    # Derivadas das funções de interpolação 
    # no centro do elemento 
    dNrs = MontadN(0.0,0.0)

    # Calcula a matriz J no ponto
    J = MontaJ(dNrs,X,Y)
    
    # Mapeia as derivadas para xy
    dNxy = CorrigedN(dNrs,J)
  
    # Monta a matriz B no ponto
    B = MontaB(dNxy)

    # Gls do elemento no sistema global
    ge = Montagls(ele,IJ)

    # Localiza o vetor de deslocamento 
    ue = U[ge]

    # A tensão σ⁰ no centro do elemento 
    σ0 = Cv*B*ue

end

#
# Rotina para calcular o carregamento adjunto e também 
# a derivada parcial 
#
function Derivada_tensao(ne,MAT,IJ,XY,U,ρ,μ,r,σY,expP,expQ)

    # Número de graus de liberdade do problema global 
    ngls = length(U)

    # Aloca o vetor de carregamento adjunto
    Fλ = zeros(ngls)

    # Aloca o vetor de derivadas parciais ∂LA / ∂ ρ_m
    ∂D = zeros(ne)

    # Matriz de Voigts
    V = SMatrix{3,3}(  [1.0 -0.5 0.0;
                       -0.5  1.0 0.0;
                        0.0  0.0 3.0] )

    
 
    # Aloca alguns arrays para não alocar no loop 
    um = Vector{Float64}(undef,8)
    gm = Vector{Int64}(undef,8)
   
    # Loop m 
    for m in LinearIndices(∂D)

        # Recupera os dados do material 
        E = MAT[m,1]
        ν = MAT[m,2]

        # Monta a matriz constitutiva
        Cv = MontaCEPT(E,ν)

        # Coordenadas dos nós do elemento 
        X,Y = MontaXY(m,IJ,XY) 

        # Derivadas das funções de interpolação 
        # no centro do elemento 
        dNrs = MontadN(0.0,0.0)

        # Calcula a matriz J no ponto
        J = MontaJ(dNrs,X,Y)
    
        # Mapeia as derivadas para xy
        dNxy = CorrigedN(dNrs,J)
  
        # Monta a matriz B no ponto
        B = MontaB(dNxy)

        # Gls do elemento no sistema global
        gm .= Montagls(m,IJ)

        # Localiza o vetor de deslocamento 
        um .= U[gm]

        # A tensão σ⁰ no centro do elemento 
        σ0 = Cv*B*um

        # Calcula a parametrização QP para o elemento 
        ρσ = fσ(ρ[m],expP,expQ)

        # Calcula a derivada da parametrização de tensão
        dρσ = dfσ(ρ[m],expP,expQ)

        # Calcula a tensão parametrizada
        σ = ρσ*σ0

        # Calcula a tensão equivalente 
        σeq = Tensao_eq(σ)

        # Calcula a restrição de tensão 
        g = gσ(σeq,σY)

        # <>
        bracket = Heaviside(μ[m]/r + g)

        # Monta o termo de ∂D[m]
        ∂D[m] = (r)*bracket*(1/(σY*σeq))*transpose(σ)*V*(dρσ*σ0)

        # Contribuição para o vetor Fλ
        fm = bracket*(1/(σY*σeq))*transpose(σ)*V*ρσ*Cv*B

        # Posiciona (e soma) nas posições globais do elemento 
        Fλ[gm] .= Fλ[gm] .+ vec(fm)

    end #m

    # Devolve os vetores 
    return -(r)*Fλ, ∂D

end