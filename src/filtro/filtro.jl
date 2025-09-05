#
# Filtro espacial de pesos lineares - Método tradicional 
# sem usar a matriz de filtragem
#
function Filtro(ne,vizinhos,pesos,a::Vector)

	# Cria um vetor de valores filtrados
	a_filt = zeros(ne)
	
	# Loop pelos elementos 
	for ele in LinearIndices(a_filt)
	
		# Vizinhos e pesos deste elemento 
		viz = vizinhos[ele]
		w   = pesos[ele]
		
		# Parte de cima da fração 
		cima = sum( a[viz].*w )
		
		# Parte de baixo da fração 
		baixo = sum(w)
	
	    # Valor filtrado 
		a_filt[ele] = cima/baixo
	
	end #ele 

    # Devolve o vetor filtrado 
	return a_filt

end


#
# Monta a matriz de filtragem
#
function Matriz_F(ne,vizinhos,pesos)

	# Aloca a matriz de filtragem
	F = spzeros(ne,ne)

	# Loop pelos elementos da malha, pegando os 
	# vizinhos e os pesos
	for ele=1:ne

		# Vizinhos do elemento
		viz = vizinhos[ele]

		# pesos de cada vizinho 
		w = pesos[ele]

		# soma dos pesos 
		soma = sum(w)

		# Coloca w/soma em cada uma das colunas 
		# associadas aos vizinhos
		F[ele,viz] .= w./soma

	end #ele

    # Retorna a matriz de filtragem
	return F

end
