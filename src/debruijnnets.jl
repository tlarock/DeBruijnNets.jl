using Graphs, SimpleWeightedGraphs, StatsBase

"""
    Convert a dictionary keyed by edge tuples to
    integer weights into 3 vectors for input
    to SimpleWeightedDigraph.

"""
function vectors_from_edgedict(edgedict::Dict{Tuple, Integer})
    src = zeros(Int64, length(edgedict))
    dest = zeros(Int64, length(edgedict))
    weights = zeros(Int64, length(edgedict))
    idx = 1
    for tup in edgedict
        src[idx] = tup[1][1]
        dest[idx] = tup[1][2]
        weights[idx] = tup[2][1]
        idx += 1
    end
    return src, dest, weights
end


"""
    Add a node to a mapping dictionary
    and return the incremented index
    for convenience.

"""
function add_to_mapping!(mapping, node, idx)
    if !haskey(mapping, node)
        mapping[node] = idx
        idx = idx + 1
    end
    return idx
end


"""
"""
function remove_selfloops!(w)
    i = 2
    e = length(w)+1
    while i < e
        if w[i-1] == w[i]
            deleteat!(w, i-1)
            e-=1
        else
            i+=1
        end
    end
end

"""
    Compute a DeBruijn graph of order k and it's
    associated first order graph from an ngram file
    with or without frequency values.

"""
function from_ngram(file::String, frequency::Bool, k::Integer) 
    # Initialize empty mappings and edgelists
    # fo = first-order
    # ko = kth-order
    fo_map = Dict{String, Integer}()
    ko_map = Dict{Tuple, Integer}()
    fo_edgelist = Dict{Tuple, Integer}()
    ko_edgelist = Dict{Tuple, Integer}()
    fin = open(file)
    fo_idx = 1
    ko_idx = 1
    for line in eachline(fin)
        s = split(line, ',')
        if frequency
            walk = s[1:end-1]
            freq = parse(Float64, s[end])
        else
            walk = s
            freq = 1
        end
        remove_selfloops!(walk) 
        for i in range(1, length(walk))
            if i > 1
                u = walk[i-1]
                v = walk[i]
		if u == v
		    println("Self-loop: $u, $v")
		end
                fo_idx = add_to_mapping!(fo_map, u, fo_idx)
                fo_idx = add_to_mapping!(fo_map, v, fo_idx)
                if haskey(fo_edgelist, (fo_map[u],fo_map[v]))
                    fo_edgelist[(fo_map[u],fo_map[v])] += freq
                else
                    fo_edgelist[(fo_map[u],fo_map[v])] = freq
                end
            end
            if i < length(walk)-k+1
                u = Tuple(walk[i:i+k-1])
                v = Tuple(walk[i+1:i+k])
		if u == v
		    println("Self-loop: $u, $v")
		end
                ko_idx = add_to_mapping!(ko_map, u, ko_idx)
                ko_idx = add_to_mapping!(ko_map, v, ko_idx)

                if haskey(ko_edgelist, (ko_map[u],ko_map[v]))
                    ko_edgelist[(ko_map[u],ko_map[v])] += freq
                else
                    ko_edgelist[(ko_map[u],ko_map[v])] = freq
                end
            end
        end
    end
    close(fin)
    
    # Convert edge dictionary representations into SimpleWeightedDigraphs
    fo_src, fo_dest, fo_weights = vectors_from_edgedict(fo_edgelist)
    forder = SimpleWeightedDiGraph(fo_src, fo_dest, fo_weights)
    ko_src, ko_dest, ko_weights = vectors_from_edgedict(ko_edgelist)
    korder = SimpleWeightedDiGraph(ko_src, ko_dest, ko_weights)
    return forder, fo_map, korder, ko_map
end


"""
Convert edges of ko to walks of length k. Use edge
weights in ko as weights for walks.
"""
function walks_from_edges(ko, ko_map, fo_map)
    A = adjacency_matrix(ko)
    cids = findall(>(0), A)
    edges_as_walks = Vector{Tuple}(undef, length(cids))
    ko_rev_map = Dict(value => key for (key, value) in ko_map)
    weights = Vector{Int64}(undef, ne(ko))
    Threads.@threads for aidx in range(1, length(cids))
        cidx = cids[aidx]
        u, v = Tuple(cidx)
        ko_u = Tuple([fo_map[j] for j in ko_rev_map[u]])
        ko_v = Tuple([fo_map[j] for j in ko_rev_map[v]])
        walk = Tuple{Vararg{Int64}}([ko_u...,ko_v[end]])
        edges_as_walks[aidx] = walk
        weights[aidx] = A[u, v]
    end
    return edges_as_walks, weights
end

"""
Construct a graph from input walks.

NOTE: Expects nodes to be integers.
"""
function from_walks(walks::Vector{Tuple}, k::Int64, rev_fo_map::Dict{Integer, String})
    ko_map = Dict{Tuple, Integer}()
    fo_edgelist = Dict{Tuple, Integer}()
    ko_edgelist = Dict{Tuple, Integer}()
    new_konode_idx = 1
    for w in walks
        for i in range(2, length(w))
            edge = w[i-1]...,w[i]
            if !haskey(fo_edgelist, edge)
                fo_edgelist[edge] = 0
            fo_edgelist[edge] += 1
            end
        end

        for i in range(1, length(w)-k)
            mapped = Tuple{Vararg{String}}([String(rev_fo_map[u]) for u in w[i:i+k]])
            u = mapped[i:i+k-1]
            v = mapped[i+1:i+k]
            if !haskey(ko_map, u)
                ko_map[u] = new_konode_idx
                new_konode_idx += 1
            end
            if !haskey(ko_map, v)
                ko_map[v] = new_konode_idx
                new_konode_idx += 1
            end
            edge = Tuple{Integer, Integer}([ko_map[u],ko_map[v]])
            if !haskey(ko_edgelist, edge)
               ko_edgelist[edge] = 0
            end
            ko_edgelist[edge] += 1
        end
    end
    fo_src, fo_dest, fo_weights = vectors_from_edgedict(fo_edgelist)
    forder = SimpleWeightedDiGraph(fo_src, fo_dest, fo_weights)
    ko_src, ko_dest, ko_weights = vectors_from_edgedict(ko_edgelist)
    korder = SimpleWeightedDiGraph(ko_src, ko_dest, ko_weights)
    return forder, korder, ko_map
end

"""
Return all k-edge walks based on first-order
graph G.
"""
function get_all_walks(G, k)
    adjacency_dict = get_adj_dict(G, adjacency_matrix(G))
    walklist = Vector{Tuple}()
    for node in range(1, nv(G))
        nodewalks = all_walks(adjacency_dict, node, k)
        append!(walklist, nodewalks)
    end
    return walklist
end

"""
Return all k-edge walks based on first-order
graph G and also return all k-order nodes.
"""
function get_all_walks_with_nodes(G, k)
    adjacency_dict = get_adj_dict(G, adjacency_matrix(G))
    walklist = Vector{Tuple}()
    nodes = Set{Tuple}()
    for node in range(1, nv(G))
        nodewalks = all_walks(adjacency_dict, node, k)
        append!(walklist, nodewalks)
        for walk in nodewalks
            u = walk[1:end-1]
            v = walk[2:end]
            push!(nodes, u)
            push!(nodes, v)
        end
    end
    return walklist, nodes
end
