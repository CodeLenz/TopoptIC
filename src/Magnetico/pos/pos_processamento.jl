#
# Calcula o campo vetorial H baseado na solução nodal 
# Φ
#
function Calcula_HB(ne,coord,connect,materials,Φ)

   # Aloca matrizes de saída em que 
   # cada linha é um elemento, a primeira 
   # coluna é ∂x e a segunda é ∂y
   H = zeros(ne,2)
   B = zeros(ne,2)

   # Loop pelos elementos, extraindo o ϕ_e
   # e montando a matriz Be no centro do elemento
   for ele=1:ne

       # Material
       mat = connect[ele,2] 

       # Find material μ
       μ = materials[mat,1]

       # Descobre o tipo de elemento 
       etype = connect[ele,1]   

       # Descobre os nós e coordenadas do elemento 
       nos, X = Nos_Coordenadas(ele,etype,coord,connect)

       # Pega somente os valores nos nós do elemento 
       ϕe = Φ[nos]

       # Monta a matriz Be no centro do elemento 
       Be,_ = Matriz_B_bi4(0.0,0.0,X)

       # Calcula Hele 
       Hele = -Be*ϕe

       # Calcula Bele
       Bele = μ*Hele

       # Armazena nas linhas das matrizes globais
       H[ele,:] .= Hele[:]
       B[ele,:] .= Bele[:]

   end

   # Retorna as matrizes com os campos vetoriais 
   # pós-processados
   return H, B

end