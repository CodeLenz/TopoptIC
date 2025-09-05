#
# Minimização de Compliance (flexibilidade) com restrição 
# de volume utilizando o critério de ótimo, filtro espacial 
# e projeção Heaviside suavizada.
#
function MinComplianceLA(arquivo,R=0.15)

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

    # Monta o vetor com as áreas de cada elemento da malha
    A = [Area_elemento(e,IJ,XY) for e=1:ne]
    
    # Com isso, o valor limite será 
    Vlimite = vf*sum(A)

    # Visualização dos resultados
    Lgmsh_export_init("saida.pos",nn,ne,XY,etypes,IJ) 

    ###########################################################
    # ------------- COMEÇO DA MODIFICAÇÃO DO LA ---------------
    ###########################################################

    # Definições para o LA
    μ = zeros(1)
    r = 1.0
 
    # Vetor de normalização do objetivo
    VALS = [0.0]

    # Número de iterações externas no LA
    niter_LA = 1


    # Podemos definir os drivers
    LA(x) = Driver_compliance(x,r,μ,MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,exP,
                              nf,FC,np,P,na,AP,A,Vlimite,VALS,"LA")

    dLA(x) = Driver_compliance(x,r,μ,MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,exP,
                               nf,FC,np,P,na,AP,A,Vlimite,VALS,"dLA")
                          
    g(x) = Driver_compliance(x,r,μ,MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,exP,
                               nf,FC,np,P,na,AP,A,Vlimite,VALS,"g")
   
    # Restrições laterais
    ci = zeros(ne)
    cs = ones(ne)                           

    #
    #                   Loop externo do LA
    #
    for k=1:niter_LA
        
        println("Iteração ", k)

        # Chama o otimizador interno
        options = WallE.Init()
        options["NITER"] = 1000
        output = WallE.Solve(LA,dLA,x,ci,cs,options)

        # Recupera a solução
        xn .= output["RESULT"]
        flag_converged = output["CONVERGED"]
        opt_norm = output["NORM"]

        # Variação das variáveis de projeto 
        dx = norm(xn.-x,Inf)
        println("||dx||   ",dx)
        println()

        # Copia a nova distribuição para a configuração atual
        x .= xn

        # Mapeamento entre x e ρ
        ρ = Mapeamento(x,MF,β,η,ρ_min)

        # Visualização das densidades
        Lgmsh_export_element_scalar("saida.pos",ρ,"ρ $k")

        # Atualiza a penalização 
        r = r*1.1  

        # Calcula as restrições no novo ponto
        gs = g(x)

        # Atualiza a estimativa de μ
        μ = r.*Heaviside.(μ/r .+ gs)

        println("Resumo")
        @show gs
        @show r
        @show μ


    end # iter
   
end