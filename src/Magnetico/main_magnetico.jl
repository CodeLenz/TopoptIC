using LinearAlgebra
using SparseArrays
using StaticArrays
using DelimitedFiles
using Lgmsh

# Carrega as rotinas do programa
include("global/auxiliar.jl")
include("elementos/bi4.jl")
include("elementos/tri.jl")
include("elementos/tetra.jl")
include("elementos/pyramid.jl")
include("elementos/hexa.jl")

include("global/carregamento.jl")
include("global/global.jl")
include("parse/parse.jl")
include("pos/pos_processamento.jl")

#
# Calcula o campo escalar Φ e também 
# os campos vetoriais H e B
#
function Analise(meshfile::String; output=true)
    
    # Le dados da malha
    nn, coord, ne, connect, materials, φm, vetor_hn, vetor_ρm, _ = Parsemsh(meshfile)

    # Precisamos de um material
    if isempty(materials)
        error("Analise:: ao menos um material é necessário")
    end

    # Calcula a matriz de rigidez global
    K = Monta_K(ne,coord,connect,materials)

    # Posições que não precisam ser calculadas no sistema de equações 
    livres = setdiff(collect(1:nn),φm)
    
    # Vetor de forças devido a h normal ao contorno
    Phn = Vetor_Phn(nn,vetor_hn,coord,connect)

    # Vetor de forças de corpo devido a ρm
    Pρm = Vetor_Pρm(nn,vetor_ρm,coord,connect)

    # Força total 
    P = Phn - Pρm

    # Aloca um vetor Φ com a dimensão completa
    Φ = zeros(nn)

    # Soluciona o sistema de equações para os gls livres
    Φ[livres] .= -K[livres,livres]\P[livres]

    # Agora podemos pós-processar a solução, calculando o H e também 
    # o B no centro de cada elemento finito
    H, B =  Calcula_HB(ne,coord,connect,materials,Φ)

    # Abre um arquivo para escrita no gmsh
    if output

        # Abre o arquivo de saída e coloca as informações
        # básicas
        etypes = connect[:,1]
        conectividades = connect[:,3:end]
        Lgmsh_export_init("saida.pos",nn,ne,coord,etypes,conectividades)

        # Grava o campo Φ
        Lgmsh_export_nodal_scalar("saida.pos",Φ,"Campo escalar") 
        
        # Grava o campo h
        Lgmsh_export_element_scalar("saida.pos",H[:,1],"Hx") 
        Lgmsh_export_element_scalar("saida.pos",H[:,2],"Hy") 

    end

    # Retorna os resultados
    return Φ, H, B
   
end