#
# Determina os vizinhos de cada elemento da malha
#
function Vizinhanca(ne,centroides,R)

    # Vetor de vetores com os vizinhos de cada elemento
	vizinhos = Vector{Vector{Int}}(undef,ne)
	
	# Vetor de vetores com os pesos de cada vizinho
	pesos = Vector{Vector{Float64}}(undef,ne)
	
	# Mínima distância entre centroides da malha
	d_min = maxintfloat(1.0) 

    # Menor e maior número de vizinhos na malha
	min_viz = ne
	max_viz = 0

	# Loop pelos elementos da malha
	for ele=1:ne
	
	    # Posição (x,y) do centróide deste elemento
		c_ele = centroides[ele,:]
		
		# Cria um vetor local para armazenar os vizinhos deste elemento 
        vetor_local = Int[]		
		
		# Cria um vetor local para armazenar os pesos dos vizinhos
        pesos_local = Float64[]		
				
		# Loop pelos elementos da malha
		for viz=1:ne
		
			# Centróide deste elementos
			c_viz = centroides[viz,:]
			
			# Distância entre os centróides
			dist = norm(c_viz.-c_ele)
			
			# Para controle 
			if viz!=ele
				d_min = min(dist,d_min)
			end

			# Se a distância for <= R, é vizinho 
			if dist<=R
			
			   # Guarda o vizinho 
			   push!(vetor_local,viz)
			   
			   # Calcula e guarda o peso 
               w = 1 - dist/R
               push!(pesos_local,w)			   
			   
			end
		
		end #viz
		
		# Armazena o vetor de vizinhos e de pesos deste elemento nos 
		# vetores de vetores 
		vizinhos[ele] = copy(vetor_local)
		pesos[ele]    = copy(pesos_local)
		
		# Mínimo e máximo número de vizinhos
		min_viz = min(min_viz, length(vetor_local))
		max_viz = max(max_viz, length(vetor_local))

	end # ele

    println("Distância mínima ",d_min)
	println("Menor número de vizinhos ",min_viz)
	println("Maior número de vizinhos ",max_viz)

    # Retorna os vizinhos e pesos
	return vizinhos, pesos

end