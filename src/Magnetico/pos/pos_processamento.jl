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
       et = connect[ele,1]   

       # Descobre os nós e coordenadas do elemento 
       nos, X = Nos_Coordenadas(ele,et,coord,connect)

       # Pega somente os valores nos nós do elemento 
       ϕe = Φ[nos]

       # Monta a matriz Be no centro do elemento 
       if et==3
          Be,_ = Matriz_B_bi4(0.0,0.0,X)
       elseif et==2
          Be = Matriz_B_tri3(X)
       elseif et==4
          Be = Matriz_B_tet4(X)
       elseif et==5
          Be,_ = Matriz_B_hex8(0.0,0.0,0.0,X)
       elseif et==7
          Be,_ = Matriz_B_pyr5(0.0,0.0,0.0,X)
       else
        error("Calcula_HB!:: Tipo de elemento não definido")
       end
       
       
       # Calcula Hele 
       Hele = -Be*ϕe

       # Calcula Bele
       # TODO Arrumar o cálculo do B (considerando o ρm)
       Bele = μ*Hele

       # Armazena nas linhas das matrizes globais
       H[ele,:] .= Hele[:]
       B[ele,:] .= Bele[:]

   end

   # Retorna as matrizes com os campos vetoriais 
   # pós-processados
   return H, B

end