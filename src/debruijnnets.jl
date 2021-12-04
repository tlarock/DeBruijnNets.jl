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
        
        for i in range(1, length(walk))
            if i > 1
                u = walk[i-1]
                v = walk[i]
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
