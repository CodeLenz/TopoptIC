

#
# Driver do problema de min volume co 
# restrição local de tensão
#
# LA(x) = Driver_tensao(x,....,"LA")
#
#
#
function Driver_tensao(x::Vector,r::Float64,μ::Vector,
                       MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,expP,expQ,
                       nf,FC,np,P,na,AP,A,σY,VALS::Vector,opcao::String)

    # Verifica a consistência de opcao
    opcao in ["volume","g","LA","dLA"] || error("Driver_tensao:: opção $opcao não é válida")

    # Mapeamento entre x e ρ
    ρ = Mapeamento(x,MF,β,η,ρ_min)

    # Matriz de rigidez global
    K = RigidezGlobal(nn,ne,MAT,ESP,XY,IJ,ρ,expP)

    # Aplica CCH na rigidez
    AplicaCCH_K!(na,AP,K) 

    # Vetor de forças concentradas global
    Fc = ForcaCglobal(nn,nf,FC)

    # Vetor de forças distribuídas global
    Ft = ForcatGlobal(nn,ESP,XY,IJ,np,P)

    # Força total 
    F = Fc .+ Ft

    # Aplica CCH na força
    AplicaCCH_F!(na,AP,F) 

    # Solução do sistema linear de Equações utilizando 
    # o método Cholesky (pois a matriz é simétrica e posdef)
    CK = cholesky(K)
    U = CK\F

    # Calcula o volume atual da estrutura V(ρ(x))
    volume = dot(ρ,A)

    # Se o usuário quer só ver o valor do volume 
    if opcao=="volume"
        return volume
    end

   # @show volume, norm(ρ)

    # Calcula as restrições de tensão 
    g = zeros(ne)

    # Loop para calcular a restrição de tensão em cada elemento 
    for ele in LinearIndices(g)

        # Tensão no centro do elemento 
        σ0 = Tensao0_centro_elemento(ele,MAT,IJ,XY,U)

        # Parametriza a tensão 
        ρσ = fσ(ρ[ele],expP,expQ)

        # Calcula a tensão equivalente no ponto 
        σeq = Tensao_eq(ρσ*σ0)

        # Calcula a restrição e armazena em g[ele]
        g[ele] = gσ(σeq,σY)

    end

    # Se o usuário quer só ver o vetor de restrições
    if opcao=="g"
        return g
    end

    # Normalização 
    if VALS[1]==0.0

        # Primeira vez que chamamos o Driver
        VALS[1] = volume
    end

    # Com essas informações, podemos calcular o LA
    LA1 = volume/VALS[1] 
    LA2 = (r/(2))*sum(Heaviside.(μ/r .+ g).^2)

    #@show LA2 

    LA = LA1 + LA2

    # Se o usuário quer o valor do LA 
    if opcao=="LA"
        return LA
    end

    #
    # Monta o vetor de "carregamento" adjunto e já aproveita
    # para calcular a derivada parcial 
    #
    Fλ, ∂D = Derivada_tensao(ne,MAT,IJ,XY,U,ρ,μ,r,σY,expP,expQ)
    
    # Aplica as cc em Fλ
    AplicaCCH_F!(na,AP,Fλ) 

    # Resolve o problema adjunto
    λ = CK\Fλ

    # Agora podemos calcular a parcela da derivada que depende de λ
    dλ = Derivada_termo_λ(ne,ρ,expP,U,λ,MAT,ESP,XY,IJ)

    # Derivada do "volume" em relação a ρ
    dV = A 

    # Assim, a derivada do LA em relação a ρ
    dLAρ = dV/VALS[1] .+ ∂D .+ dλ

    # Corrige a derivada para ser em relação a x
    return dMapeamento(dLAρ,MF,x,β,η,ρ_min)

end