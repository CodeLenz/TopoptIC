#
# Read a .msh file and process data to a specific format
#
#
# Physical groups que esta rotina consegue processar
#
#
# Material,nome,id,μ [ surfaces (and volumes) ]
#
# Open [ lines and/or nodes and/or surfaces]
# 
# bn,value [ lines and/or nodes and/or surfaces]
#
#
# Elementos que estão implementados
#
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
    materials = Dict{String,Union{Float64,Int64,Vector{Int64}}}[]

    # Local dict inside the loop
    localD_m = Dict{String,Union{Float64,Int64,Vector{Int64}}}()

    # Vector with Dicts of normal b on nodes
    bn = Dict{String,Union{Float64,Matrix{Int64}}}[]

    # Local dict inside the loop
    localD_bn = Dict{String,Union{Float64,Matrix{Int64}}}()

    # Vector of OPEN nodes
    nodes_open = Int64[]

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
            push!(materials,copy(localD_m))

      elseif  occursin("Open",st[1])

            # Find nodes 
            nodes = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            nodes_open = vcat(nodes_open,nodes)

      
      elseif  occursin("bn",st[1])

            # Clean dictionary to store local data
            empty!(localD_bn)

            # Normal b
            localD_bn["value"] = parse(Float64,st[2])

            # Find nodes 
            nodes_bn = Lgmsh.Readnodesgroup(meshfile,name)

            # If 2D  - Find element and edges
            # else   - Find element faces
            # Vamos continuar chamando de edges, mesmo em 3D
            eleedges = Int64[]
            edges = Int64[]
            for tt in et
                if dimensao==2
                   eleedges_,edges_ = FindElementsEdges(tt,ne,etypes,connect,nodes_bn)
                else
                    eleedges_,edges_ = FindElementsFaces(tt,ne,etypes,connect,nodes_bn)
                end
                if !isempty(eleedges_)
                    push!(eleedges,eleedges_...)
                    push!(edges,edges_...)
                end
            end

            # Append
            localD_bn["elements"] = [eleedges edges]

            # Copy the dict to the vector of velocities
            push!(bn,copy(localD_bn))

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
    for mat in materials

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

    # Testing...
    #=
    Mesh(nn, coord, ne, connect2, materials2, unique!(nodes_open), velocities, unique!(nodes_pressure), 
    pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed)
    =#

    # Return processed data
    return nn, coord, ne, connect2, materials2, unique!(nodes_open), bn, centroids

end    