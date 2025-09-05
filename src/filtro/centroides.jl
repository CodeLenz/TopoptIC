#
# Monta uma matriz com os centróides dos elementos da malha
#
function Centroides(ne,IJ,XY)

	# Inicializa a matriz de centróides ne x 2
	centroides = zeros(ne,2)
	
	# Loop pelos elementos
	for ele=1:ne
	
		# Nós do elemento 
		nos = IJ[ele,:]
		
		# Coordenadas X dos nós do elemento
		x = XY[nos,1]
	
	    # Média dos valores em x
        xm = mean(x)		
		
	    # Coordenadas Y dos nós do elemento
		y = XY[nos,2]
	
	    # Média dos valores em y
        ym = mean(y)		
		
		# Grava na linha da matriz de centróides
		centroides[ele,1] = xm 
		centroides[ele,2] = ym 
		
		    
	end # ele

    # Retorna a matriz de centróides
	return centroides

end