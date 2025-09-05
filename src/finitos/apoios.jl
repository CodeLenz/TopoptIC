#
# Aplica as condições de contorno essenciais
# Homogêneas
#
function AplicaCCH_K!(na,AP,K)

    # Laço sobre as linhas de AP
    for l = 1:na

        # Nó
        no = AP[l,1]

        # gll
        gll = AP[l,2]

        # glg
        glg = Int(2*(no-1)+gll)

        # Zera linha
        K[glg,:] .= 0

        # Zera coluna
        K[:,glg] .= 0

        # Coloca 1 na diagonal da matriz
        K[glg,glg] = 1

    end #l

end


#
# Aplica as condições de contorno essenciais
# Homogêneas
#
function AplicaCCH_F!(na,AP,F)

    # Laço sobre as linhas de AP
    for l = 1:na

        # Nó
        no = AP[l,1]

        # gll
        gll = AP[l,2]

        # glg
        glg = Int(2*(no-1)+gll)

        # Zera a linha em F
        F[glg] = 0

    end #l

end