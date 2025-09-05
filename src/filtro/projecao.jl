#
# Projeção Heaviside
#
function Projecao(a,β,η)
  
    # Cria vetor de saída
    b = similar(a)

    # Ctes fora do loop 
    c1 = tanh(β*η)
    c2 = tanh(β*(1-η))
    c12 = c1+c2

    # Para cada posição de a, gera a Projeção b
    for i in LinearIndices(a)

        c3 = tanh(β*(a[i]-η)) 

        b[i] = (c1+c3)/c12

    end 

    # Retorna as variáveis projetadas
    return b

end


#
# Correção da derivada da projeção
#
function dProjecao(a,β,η)
  
    # Cria vetor de saída
    db = similar(a)

    # Cts fora do loop
    c1 = tanh(β*η)
    c2 = tanh(β*(1-η))
   
    # Para cada posição de a, gera a derivada da 
    # projeção b
    for i in LinearIndices(a)

        c3 = sech(β*(a[i]-η)) 

        db[i] = (β*c3^2)/(c1+c2)

    end 

    # Retorna as correções de derivadas devido a projeção
    return Diagonal(db)
    
end