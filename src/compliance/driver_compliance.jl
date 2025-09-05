
#
# Driver do problema de min compliance com 
# restrição de volume
#
# LA(x) = Driver_compliance(x,....,"LA")
#
#
#
function Driver_compliance(x::Vector,r::Float64,μ::Vector,
                           MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,exP,
                           nf,FC,np,P,na,AP,A,Vlimite,VALS::Vector,opcao::String)

    # Verifica a consistência de opcao
    opcao in ["compliance","g","LA","dLA"] || error("Driver_compliance:: opção $opcao não é válida")

    # Mapeamento entre x e ρ
    ρ = Mapeamento(x,MF,β,η,ρ_min)

    # Matriz de rigidez global
    K = RigidezGlobal(nn,ne,MAT,ESP,XY,IJ,ρ,exP)

    # Aplica CCH na rigidez
    AplicaCCH_K!(na,AP,K) 

    # Vetor de forças concentradas global
    Fc = ForcaCglobal(nn,nf,FC)

    # Vetor de forças distribuídas global
    Fb = ForcatGlobal(nn,ESP,XY,IJ,np,P)

    # Força total 
    F = Fc .+ Fb

    # Aplica CCH na força
    AplicaCCH_F!(na,AP,F) 

    # Solução do sistema linear de Equações utilizando 
    # o método Cholesky (pois a matriz é simétrica e posdef)
    CK = cholesky(K)
    U = CK\F

    # Calcula a compliance
    compliance = transpose(F)*U

    # Se o usuário quer só ver o valor da compliance
    if opcao=="compliance"
        return compliance
    end

    # Calcula o volume atual da estrutura V(ρ(x))
    volume = dot(ρ,A)

    # Só temos uma restrição 
    g = [volume/Vlimite - 1]

    # Se o usuário quer só ver o vetor de restrições
    if opcao=="g"
        return g
    end

    # Normalização 
    if VALS[1]==0.0
        # Primeira vez que chamamos o Driver
        VALS[1] = compliance
    end

    # Com essas informações, podemos calcular o LA
    LA1 = compliance/VALS[1] 
    LA2 = (r/2)*Heaviside(μ[1]/r + g[1])^2
    LA = LA1 + LA2

    # if opcao=="dLA"
    #    @show LA1, LA2
    # end

    # Se o usuário quer o valor do LA 
    if opcao=="LA"
        return LA
    end

    #
    #
    #
    # Monta o vetor de "carregamento" adjunto
    #
    Fλ = -F/(VALS[1]) # <- depende do problema

    # Aplica as cc em Fλ
    AplicaCCH_F!(na,AP,Fλ) 

    # Resolve o problema adjunto
    λ = CK\Fλ

    # Agora podemos calcular a parcela da derivada que depende de λ
    dλ = Derivada_termo_λ(ne,ρ,exP,U,λ,MAT,ESP,XY,IJ)

    # Derivada do "volume" em relação a ρ
    dV = A 

    # Com isso, podemos montar a derivada do LA em relação a ρ
    dLAρ =  r*Heaviside(μ[1]/r + g[1]).*(dV/Vlimite) .+ dλ

    # Corrige a derivada para ser em relação a x
    return dMapeamento(dLAρ,MF,x,β,η,ρ_min)

end