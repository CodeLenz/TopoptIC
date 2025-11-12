#
# Rotina que monta o vetor de "forças" de contorno 
# 
#
function Vetor_Phn(nn,hn,coord,connect)

  # Aloca o vetor de carregamento
  P = zeros(nn)

  # Loop pelo vetor bn. Cada linha deste vetor é um dicionário
  for dhn in hn

    # Recover data from Dictionary
    valor = dhn["value"]

    # Elementos em que a condição de contorno está sendo imposta
    elements   = dhn["elements"] 

    # Agora precisamos fazer um loop sobre os elementos
    # e suas arestas
    for i in axes(elements,1)

      # Element and edge
      ele  = elements[i,1]
      edge = elements[i,2]

      # Element type
      et = connect[ele,1]

      # Find nodes and coordinates
      nos,X = Nos_Coordenadas(ele,et,coord,connect)

      # value
      val = valor

      # Local vector 
      #Pn = Edge_load_local_bi4(edge,val,X)

      if et==3
        Pn = Edge_load_local_bi4(edge,val,X)
      elseif et==2
        Pn = Edge_load_local_tri3(edge,val,X)
      elseif et==4
        Pn = Face_load_local_tet4(edge,val,X)
      elseif et==5
        Pn = Face_load_local_hex8(edge,val,X)
      elseif et==7
        Pn = Face_load_local_pyr5(edge,val,X)
      else
        error("Vetor_P!:: Tipo de elemento não definido")
      end
      

      # Add to the global vector
      P[nos] .+= Pn

    end #i
 
  end # dict

  # Retorna o vetor de carregamentos
  return P

end

#
#
# Rotina que monta o vetor de "forças" de corpo
# 
#
function Vetor_Pρm(nn,vetor_ρm,coord,connect)

  # Aloca o vetor de carregamento
  P = zeros(nn)

  # Loop pelo vetor bn. Cada linha deste vetor é um dicionário
  for dρm in vetor_ρm

    # Recover data from Dictionary
    valor = dρm["value"]

    # Elementos em que a condição de contorno está sendo imposta
    elements   = dρm["elements"] 

    # Agora precisamos fazer um loop sobre os elementos
    # e suas arestas
    for i in axes(elements,1)

      # Element 
      ele  = elements[i]

      # Element type
      et = connect[ele,1]

      # Find nodes and coordinates
      nos,X = Nos_Coordenadas(ele,et,coord,connect)

      # Local vector
      if et==3
        Pb = Body_load_local_bi4(valor,X)
      elseif et==2
        Pb = Body_load_local_tri3(valor,X)
      elseif et==4
        Pb = Body_load_local_tet4(valor,X)
      elseif et==5
        Pb = Body_load_local_hex8(valor,X)
      elseif et==7
        Pb = Body_load_local_pyr5(valor,X)
      else
        error("Vetor_P!:: Tipo de elemento não definido")
      end
      

      # Add to the global vector
      P[nos] .+= Pb

    end #i
 
  end # dict

  # Retorna o vetor de carregamentos
  return P

end