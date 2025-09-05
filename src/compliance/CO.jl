#
# Critério de ótimo para Compliance
#
function CO(ne,dC,dV,x,Vlimite,A,MF,β,η,ρ_min,δ=0.1,tol=1E-6)

    # Inicializa o intervalor do
    # multiplicador de Lagrange
    λl = 0.0
    λm = maxintfloat(1.0)

    # Inicializa o vetor de variáveis de projeto
    # atualizadas
    xn = copy(x)

    # Loop da biseção
    while(λm-λl>tol)

        # Biseção do λ
        λ = (λl + λm)/2
        
        # Loop para ajustar o β
        for m=1:ne

            # Calcula o βρ
            βρ = - dC[m]/( λ*dV[m] )

            # Escalonamento do ρ
            can = x[m]*sqrt(βρ)

            # Testamos pelo limite móvel
            cansup = min(x[m]+δ,1)
            caninf = max(x[m]-δ,0)

            # Satisfaz os limites móveis
            # e as restrições laterais
            xn[m] = max(min(can,cansup),caninf)

        end #m

        # Calcula o volume considerando a nova distribuição de 
        # material. Vamos ter que fazer uma projeção aqui...
        ρ = Mapeamento(xn,MF,β,η,ρ_min)
        Vnovo = Volume(ρ,A)

        # Lógica da biseção
        if Vnovo-Vlimite > 0
            λl = λ
        else
            λm = λ
        end

    end # biseção

    # Retorna a nova distribuição das variáveis de projeto
    return xn

end