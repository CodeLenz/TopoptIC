

#
# Minimização de volume com restrição de tensão local, filtro espacial 
# e projeção Heaviside suavizada.
#
function MinVolσ(arquivo,R=0.15; verifica_derivada=false)


    # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca
    # do gmsh
    if occursin(".geo",arquivo)
       
       # Gera a malha
       gmsh.initialize()
       gmsh.open(arquivo)
       gmsh.model.mesh.generate(2)
       
       # Cria o mesmo nome, mas com .msh
       mshfile = replace(arquivo,".geo"=>".msh")

       # Cria o .msh
       gmsh.write(mshfile)
      
    else 

       # Assumimos que já passaram o .msh (seria bom testar...)
       mshfile = arquivo

    end

    # Nome do arquivo .pos (sem o diretório)
    # Recupera o nome do arquivo sem caminho
    nome_msh = basename(mshfile)
    nomepos  = replace(nome_msh,".msh"=>".pos")


    # Entrada de dados
    nn,XY,ne,IJ,MAT,ESP,nf,FC,np,P,na,AP,nfb,FB,etypes,centroides = ConversorFEM1(mshfile)

    # Determina a vizinhança de cada elemento da malha
	vizinhos,pesos = Vizinhanca(ne,centroides,R)

    # Monta a matriz de filtragem 
    MF =  Matriz_F(ne,vizinhos,pesos)

    #
    #                      Processamento
    #

    # Define a distribuição inicial de variáveis de projeto
    # Variáveis matemáticas
    x = ones(ne)

    # Define um vetor de atualização 
    xn = similar(x)

    # Define o expoente da parametrização SIMP
    expP = 3.0

    # Define o expoente da parametrização QP
    expQ = 1.5

    # Tensão de escoamento 
    σY = 30E6

    # Fator de ajuste da projeção Heaviside
    β = 5.0

    # Centro do Heaviside
    η = 0.5

    # Densidade relativa mínima
    ρ_min = 1E-3

    # Monta o vetor com as áreas de cada elemento da malha
    A = [Area_elemento(e,IJ,XY) for e=1:ne]
    
    # Visualização dos resultados
    Lgmsh_export_init(nomepos,nn,ne,XY,etypes,IJ) 

    ###########################################################
    # ------------- COMEÇO DA MODIFICAÇÃO DO LA ---------------
    ###########################################################

    # Definições para o LA
    μ = zeros(ne)
    r = 10.0
 
    # Vetor de normalização do objetivo
    VALS = [0.0]

    # Número de iterações externas no LA
    niter_LA = 10

    # Podemos definir os drivers
    LA(x) = Driver_tensao(x,r,μ,MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,expP,expQ,
                              nf,FC,np,P,na,AP,A,σY,VALS,"LA")

    dLA(x) = Driver_tensao(x,r,μ,MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,expP,expQ,
                              nf,FC,np,P,na,AP,A,σY,VALS,"dLA")
                          
    g(x) = Driver_tensao(x,r,μ,MF,β,η,ρ_min,nn,ne,MAT,ESP,XY,IJ,expP,expQ,
                              nf,FC,np,P,na,AP,A,σY,VALS,"g")
   

    # Restrições laterais
    ci = zeros(ne)
    cs = ones(ne)                           


    #
    # Verificação de derivada
    #
    if verifica_derivada

        # Desliga a normalização do objetivo no LA
        VALS = [1.0]

        # Vamos testar a derivada
        xa = rand(ne)

        # Calcula as restrições no ponto atual 
        gxa = g(xa)

        # Calcula a derivada "analítica"
        dfr = dLA(xa)

        # Calcula a derivada numérica
        dfn = DiferencasFinitas(LA,xa,1E-6)

        return xa, gxa, dfr, dfn

    end

    #
    #                   Loop externo do LA
    #
    for k=1:niter_LA
        
        println("Iteração ", k)

        # Chama o otimizador interno
        options = WallE.Init()
        options["NITER"] = 100
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
        Lgmsh_export_element_scalar(nomepos,ρ,"ρ $k")

        # Atualiza a penalização 
        r = r*1.1  

        # Calcula as restrições no novo ponto
        gs = g(x)

        # Atualiza a estimativa de μ
        μ = r.*Heaviside.(μ/r .+ gs)

        println("Resumo")
        @show gs[gs.>0]
        @show r
        @show μ[μ.>0]


    end # iter
   
end