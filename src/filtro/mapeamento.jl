#
# Mapeamento entre x e ρ
#
function Mapeamento(x,MF,β,η,ρ_min)

    # Primeira operação é a filtragem 
    a = MF*x

    # Depois a projeção
    b = Projecao(a,β,η)

    # Terceira operação é o offset 
    ρ = ρ_min .+ (1-ρ_min).*b

end

#
# Correção do gradiente de ρ -> x
#
# df -> derivada em relação a ρ
# 
#
function dMapeamento(df,MF,x,β,η,ρ_min)

    # Faz o filtro x->a 
    a = MF*x

    # Cria a matriz R (matriz diagonal)
    R = dProjecao(a,β,η)

    # Faz a correção
    (1-ρ_min).*R*transpose(MF)*df

end

