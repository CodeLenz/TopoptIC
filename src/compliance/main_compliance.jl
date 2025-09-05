#
# Minimização de Compliance (flexibilidade) com restrição 
# de volume utilizando o critério de ótimo, filtro espacial 
# e projeção Heaviside suavizada.
#
function MinCompliance(arquivo,R=0.15)

    # Entrada de dados
    nn,XY,ne,IJ,MAT,ESP,nf,FC,np,P,na,AP,nfb,FB,etypes = ConversorFEM1(arquivo)

    # Calcula os centróides dos elementos 
    centroides = Centroides(ne,IJ,XY)	

    # Determina a vizinhança de cada elemento da malha
	vizinhos,pesos = Vizinhanca(ne,centroides,R)

    # Monta a matriz de filtragem 
    MF =  Matriz_F(ne,vizinhos,pesos)

    #
    #                      Processamento
    #

    # Fração de volume
    vf = 0.3

    # Define a distribuição inicial de variáveis de projeto
    # Variáveis matemáticas
    x = vf*ones(ne)

    # Define um vetor de atualização 
    xn = similar(x)

    # Define o expoente da parametrização SIMP
    exP = 3.0

    # Fator de ajuste da projeção Heaviside
    β = 1.0

    # Centro do Heaviside
    η = 0.5

    # Densidade relativa mínima
    ρ_min = 1E-3

    # Número de iterações do loop de otimização sequencial
    niter = 100

    # Cria as variáveis que queremos ver depois do loop
    U  = zeros(2*nn)
    dC = zeros(ne)
    dV = zeros(ne)

    # Monta o vetor com as áreas de cada elemento da malha
    A = [Area_elemento(e,IJ,XY) for e=1:ne]
    
    # Com isso, o valor limite será 
    Vlimite = vf*sum(A)

    # Visualização dos resultados
    Lgmsh_export_init("saida.pos",nn,ne,XY,etypes,IJ ) 

    #
    #                   Loop da otimização sequencial
    #
    for iter=1:niter

        # Mapeamento entre x e ρ
        ρ = Mapeamento(x,MF,β,η,ρ_min)

        # Visualização das densidades
        Lgmsh_export_element_scalar("saida.pos",ρ,"ρ $iter")

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
        U .= CK\F

        # Calcula a compliance
        compliance = transpose(F)*U 

        # Sumário da iteração 
        println("Iteração ", iter) 
        println("Objetivo ", compliance)
        println("Volume   ", dot(ρ,A), " com limite ", Vlimite )

        #Lgmsh_export_nodal_vector("saida.pos",U,2,"Deslocamentos")

        # Calcula a derivada da flexibilidade em relação a ρ
        dC .= Derivada_C(ne,ρ,exP,U,MAT,ESP,XY,IJ)

        # Derivada do "volume" em relação a ρ
        dV .= A

        # Corrige as derivadas devido ao mapeamento x->ρ
        dC_x = dMapeamento(dC,MF,x,β,η,ρ_min)
        dV_x = dMapeamento(dV,MF,x,β,η,ρ_min)

        # Atualiza as variáveis de projeto usando o CO
        xn .= CO(ne,dC_x,dV_x,x,Vlimite,A,MF,β,η,ρ_min)

        # Variação das variáveis de projeto 
        dx = norm(xn.-x,Inf)
        println("||dx||   ",dx)
        println()

        # Copia a nova distribuição para a configuração atual
        x .= xn

    end # iter
   
   
    # Exporta os deslocamentos nodais
    #Lgmsh_export_nodal_vector("saida.pos",U,2,"Deslocamentos")

end