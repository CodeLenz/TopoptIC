module TopoptIC
    #
    # Programa principal
    #
    using LinearAlgebra
    using SparseArrays
    using Statistics
    using StaticArrays
    using Lgmsh
    using Gmsh
    using WallE
    using ProgressMeter

    #
    # Carrega as rotinas
    #
    include("heaviside.jl")
    
    # Derivada do λ
    include("adjunto/adjunto.jl") 

    # finitos
    include("finitos/material.jl")
    include("finitos/elemento.jl")
    include("finitos/contorno.jl")
    include("finitos/global.jl")
    include("finitos/apoios.jl")
    include("finitos/conversor.jl")
    
    # Filtragem e projeção
    include("filtro/vizinhanca.jl")
    include("filtro/filtro.jl")
    include("filtro/mapeamento.jl")
    include("filtro/projecao.jl")
    
    # Volume
    include("volume/volume.jl")

    # Tensão
    include("tensao/tensao.jl")
    include("tensao/driver_tensao.jl")
    include("tensao/main_LA_tensao.jl")

    # Diferenças finitas
    include("df/df.jl")

    # Exporta a rotina principal de análise 
    export  MinVolσ

end
