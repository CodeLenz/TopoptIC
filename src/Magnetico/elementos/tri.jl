 
#
# Funções de interpolação do triângulo
#



#
# Calcula as matrizes Ke e Me para um elemento 
#
function Ke_tri3(μ,X)

  # Mapeamento para facilitar a notação
  x1,x2,x3 = X[:,1]
  y1,y2,y3 = X[:,2]

  # Termo em comum para a rigidez
  comum = ((2*x2-2*x1)*y3+(2*x1-2*x3)*y2+(2*x3-2*x2)*y1)

  # Termos da rigidez
  k11 = (y3^2-2*y2*y3+y2^2+x3^2-2*x2*x3+x2^2)/comum
  k12 = -(y3^2+(-y2-y1)*y3+y1*y2+x3^2+(-x2-x1)*x3+x1*x2)/comum
  k13 = ((y2-y1)*y3-y2^2+y1*y2+(x2-x1)*x3-x2^2+x1*x2)/comum
  k22 = (y3^2-2*y1*y3+y1^2+x3^2-2*x1*x3+x1^2)/comum
  k23 = -((y2-y1)*y3-y1*y2+y1^2+(x2-x1)*x3-x1*x2+x1^2)/comum
  k33 = (y2^2-2*y1*y2+y1^2+x2^2-2*x1*x2+x1^2)/comum
  Ke =  μ * @SMatrix [k11 k12 k13 ; k12 k22 k23 ; k13 k23 k33]

  return Ke

end

# ===================================================================================
# Force vector for a bi3 element 
# local (normal) surface load.
#
function Edge_load_local_tri3(edge,qn,X::Matrix)

  # As we assume cte load
  # and the element is linear
  # Basic test
  edge in 1:3 || throw("Map_load_local_tri3::Invalid edge")

  F1 = 1.0
  F2 = 1.0
  F3 = 1.0

  if edge==1

    dx = X[2,1]-X[1,1]
    dy = X[2,2]-X[1,2]
    F3 = 0.0

  elseif edge==2

    # (2)->(3)
    dx = X[3,1]-X[2,1]
    dy = X[3,2]-X[2,2]
    F1 = 0.0
            
  else

    # (3)->(1)
    dx = X[3,1]-X[1,1]
    dy = X[3,2]-X[1,2]
    F2 = 0.0
    
  end

  # Comprimento da aresta
  L  = sqrt(dx^2 + dy^2)

  # Resposta
  F   = (L/2)*qn*[F1;F2;F3]
  
  # Return F
  return F

end



# ===================================================================================
# Calcula a área do elemento
#
function Area_tri3(X::Matrix)

    # Monta a matriz para o cálculo do determinante
    MA = @SMatrix    [1 X[1,1] X[1,2] ;
                      1 X[2,1] X[2,2] ;
                      1 X[3,1] X[3,2] ]

    # Retorna o determinante
    0.5*det(MA)

end


#
# Vetor de forças de corpo para o tri3
#
function Body_load_local_tri3(ρm,μ, X)

    # Aloca o vetor 3 × 1 
    Fe = zeros(3)
  
    # Integração por quadratura de Gauss-Legendre
    pg = (5/3)*ones(2)
    wg = 1/3*ones(2)
  
    @inbounds for i=1:2
        # Ponto e peso nesta dimensão
        r = pg[i]
        wr = wg[i]
  
        @inbounds for j=1:2

            # Ponto e peso nesta dimensão
            s = pg[j]
            ws = wg[j]
  
            # Calcula a função de interpolação 
            N = Matriz_N_bi4(r,s)

            # Somatórios
            Fe .= Fe + μ*ρm*N'*det(J) 
  
        end
    end
  
    # Retorna o vetor de força de corpo
    return Fe
  

end
