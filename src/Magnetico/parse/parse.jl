#
# Read a .msh file and process data to a specific format
#
#
# Physical groups que esta rotina consegue processar
#
#
# Material,nome,id,μ [ surfaces (and volumes) ]
#
# φm [ lines and/or nodes and/or surfaces]
# 
# hn,value [ lines and/or nodes and/or surfaces]
#
# ρm,value [surfaces/volumes]
#
#
#        Elementos que estão implementados
# *** ESTAMOS USANDO SÓ O TIPO 3 POR ENQUANTO ***
#          2D
# 2 -> triangular (linear)
# 3 -> quadrangular (linear)
#
#          3D
# 4 -> Tetrahedra (linear)
# 5 -> hexaedra (linear)
# 7 -> pyramid (linear)
#
#
function Parsemsh(meshfile::String,verbose=false)
 
    # Primeiro precisamos definir se a malha é 2D ou 3D
    elist = Lgmsh_import_etypes(meshfile)

    # Se tivermos elementos do 4/5/7, então é 3D. Do contrário,
    # é 2D. Observe que ter 2/3 não é uma indicação direta de 
    # que a malha é 2D, pois o gmsh também gera esses elementos
    # para malhas 3D.
    dimensao = 2
    et = [2,3]
    if (4 in elist) || (5 in elist) || (7 in elist)
        dimensao = 3
        et = [4,5,7]
    end

    if verbose
        println("Solucionando um problema de dimensão $dimensao")
    end

    # Maximum number of nodes in the elements of the mesh
    nmax = maximum(Lgmsh_nodemap()[et])

    # Read mesh
    nn, coord, ne, etypes, connect, centroids, etags = Readmesh(meshfile,et)

    # Le todos os grupos físicos do arquivo 
    pgroups, pgnames = Lgmsh_import_physical_groups(meshfile)

    # Vector with Dicts of materials
    vector_materials = Dict{String,Union{Float64,Int64,Vector{Int64}}}[]

    # Local dict inside the loop
    localD_m = Dict{String,Union{Float64,Int64,Vector{Int64}}}()

    # Vector with Dicts of materials
    vector_ρm = Dict{String,Union{Float64,Vector{Int64}}}[]

    # Local dict inside the loop
    localD_ρm = Dict{String,Union{Float64,Vector{Int64}}}()

    # Vector with Dicts of normal h on nodes
    vector_hn = Dict{String,Union{Float64,Matrix{Int64}}}[]

    # Local dict inside the loop
    localD_hn = Dict{String,Union{Float64,Matrix{Int64}}}()

    # Vector of prescribed φm nodes
    nodes_φm = Int64[]

    # Maximum id in materials
    max_id = 0

    # Loop over groups
    for g in LinearIndices(pgnames)

      # Name
      name = pgnames[g]

      # Split the string by ","
      st = split(name,",")

      # Check if Material
      if occursin("Material",st[1])

            # Clean dictionary to store local data
            empty!(localD_m)

            # We expect the following data
            # name, id, dens, c, Z
            id   = parse(Int64,st[3])
            μ    = parse(Float64,st[4])

            # Populate local dict
            localD_m["id"]   = id
            localD_m["μ"]    = μ
          
            # Store maximum id
            max_id = max(max_id,id)

            # Now we must find wich elements are associated to this group
            elems_domain = Lgmsh.Readelementsgroup(meshfile,name,etags)

            # And store in the dict in Forward order (smaller to the largest)
            localD_m["elements"] = sort(elems_domain)

            # Copy the dict to the vector of materials
            push!(vector_materials,copy(localD_m))

       # densidade de carga volumétrica
       elseif occursin("ρm",st[1])

            # Clean dictionary to store local data
            empty!(localD_ρm)

            # We expect a valune
            ρ = parse(Float64,st[2])

            # Populate local dict
            localD_ρm["value"]   = ρ
          
            # Now we must find wich elements are associated to this group
            elems_domain = Lgmsh.Readelementsgroup(meshfile,name,etags)

            # And store in the dict in Forward order (smaller to the largest)
            localD_ρm["elements"] = sort(elems_domain)

            # Copy the dict to the vector of materials
            push!(vector_ρm,copy(localD_ρm))

      elseif  occursin("φm",st[1])

            # Find nodes 
            nodes = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            nodes_φm = vcat(nodes_φm,nodes)

      
      elseif  occursin("hn",st[1])

            # Clean dictionary to store local data
            empty!(localD_hn)

            # Normal b
            localD_hn["value"] = parse(Float64,st[2])

            # Find nodes 
            nodes_hn = Lgmsh.Readnodesgroup(meshfile,name)

            # If 2D  - Find element and edges
            # else   - Find element faces
            # Vamos continuar chamando de edges, mesmo em 3D
            eleedges = Int64[]
            edges = Int64[]
            for tt in et
                if dimensao==2
                   eleedges_,edges_ = FindElementsEdges(tt,ne,etypes,connect,nodes_hn)
                else
                    eleedges_,edges_ = FindElementsFaces(tt,ne,etypes,connect,nodes_hn)
                end
                if !isempty(eleedges_)
                    push!(eleedges,eleedges_...)
                    push!(edges,edges_...)
                end
            end

            # Append
            localD_hn["elements"] = [eleedges edges]

            # Copy the dict to the vector of velocities
            push!(vector_hn,copy(localD_hn))

      end #if

    end 

    # We can now process data to build connectivities with material and
    # element types 
    connect2 = zeros(Int64,ne,nmax+2)

    # Some data are already processed
    connect2[:,1] .= etypes
    connect2[:,3:end] .= connect

    # Materials as a matrix
    materials2 = zeros(max_id,1)
   
    # loop over vector of material dicts
    for mat in vector_materials

        # id
        id = mat["id"]

        # elements
        elements = mat["elements"]
        
        # Copy
        connect2[elements,2] .= id

        v = [mat["μ"]]

        # fill line of materials2
        materials2[id,:] = v
         
    end


    # Passada para converter os elementos do tipo 2 para tipo 3
    println("CUIDADO:: convertendo elementos do tipo 2 para o tipo 3")
    ntipo2 = 0
    for i=1:ne

        # Testa se o elemento é  do tipo 2
        if connect2[i,1]==2
           
            # Acumula contador de tipo 2
            ntipo2 += 1

            # Converte para 3 
            connect2[i,1] = 3

            # Copia o nó 3 para o 4
            connect2[i,end] = connect2[i,end-1]
        end

    end

    println("Convertendo $ntipo2 de $ne elementos")

    # Return processed data
    return nn, coord, ne, connect2, materials2, sort!(unique!(nodes_φm)), vector_hn, vector_ρm, centroids

end    