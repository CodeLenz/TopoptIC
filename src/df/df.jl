#
# Calcula o gradiente de uma função f(x) no entorno do ponto 
# x0
#
function DiferencasFinitas(f::Function,x0::Vector,δ=1E-6)

    # Número de variáveis 
    n = length(x0)

    # Valor atual da função 
    f0 = f(x0)

    # Aloca o vetor gradiente 
    ∇f = zeros(n)

    # Loop pelas posições do gradiente
    @showprogress "calculando df" for i in LinearIndices(∇f)

        # Backup do valor atual de x0
        x0i = x0[i]

        # Perturba a posição i de x0 para frente
        x0[i] += δ

        # Calcula a função no novo ponto 
        ff = f(x0)

        # Armazena a derivada
        ∇f[i] = (ff-f0)/δ

        # Desfaz a perturbação
        x0[i] = x0i

    end

    # Devolve o gradiente
    return ∇f

end