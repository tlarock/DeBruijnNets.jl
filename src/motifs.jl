const alphabet = ['a','b','c','d','e','f','g',
		  'h','i','j','k','l','m','n',
		  'o','p','q','r','s','t','u',
		  'v','w','x','y','z']
"""
Project a walk into the alphabet.
"""
function alphabet_projection(walk::Tuple)
    projected = Vector{Char}()
    map = Dict()
    map_idx = 1
    for w in walk
        if !haskey(map, w)
            map[w] = alphabet[map_idx]
            map_idx += 1
        end
        push!(projected, map[w])
    end
    return Tuple(projected)
end

"""
    Accept a vector of tuples representing a set of walks. Project each walk into a motif. If no weights are provided, each has weight 1.
"""
function count_motifs(walks::Vector{Tuple}, weights=ones(length(walks)))
    motifs = Dict()
    widx = 1
    for walk in walks
        p = alphabet_projection(walk)
        if !haskey(motifs, p)
            motifs[p] = 0
        end
        motifs[p] += weights[widx]
        widx += 1
    end
    return motifs
end

"""
    Accept a vector of tuples represeting a set of walks. Project each walk into a motif and record the frequency of the motif, and which edges correspond to each motif.
"""
function motifs_with_mapping(walks::Vector{Tuple}; weights::Vector{Int64}=ones(length(walks)))
    motifs = Dict()
    motif_edges = Dict()
    widx = 1
    for walk in walks
        p = alphabet_projection(walk)
        if !haskey(motifs, p)
            motifs[p] = 0.0
            motif_edges[p] = Set()
        end

        motifs[p] += weights[widx]
        push!(motif_edges[p], widx)
        widx += 1
    end
    return motifs, motif_edges
end

"""
    Accept a vector of tuples represeting a set of walks as well as a set of orders at which to compute motifs.Walks with exactly k in k_vals edges are counted only towards
    the exact motif they map to. Walks with less than min(k_vals) edges are ignored. For walks
    with more than max(k_vals) edges we slide a max(k_vals) window.
 Project each walk into a motif and record the frequency of the motif, and which edges correspond to each motif.
"""
function exclusive_motifs(walks::Vector{Tuple}, k_vals::Vector{Int64}; weights::Vector{Int64}=ones(length(walks)))
    k_set = Set(k_vals)
    max_k = maximum(k_set)
    motifs = Dict()
    motif_edges = Dict()
    widx = 1
    for walk in walks
        if (length(walk)-1) in k_set
            # Compute only 1 motif for walks with lengths in the k_set
            walk_motifs = [(alphabet_projection(walk), walk)]
        elseif (length(walk)-1) > max_k
            # Slide window of maximum length for walks with lengths
            # larger than the maximum k
            walk_motifs = Vector()
            for i in 1:length(walk)-max_k
                push!(walk_motifs, (alphabet_projection(walk[i:i+max_k]), walk[i:i+max_k]))
            end
        else
            # Ignore walks less than the minimum length
            widx += 1
            continue
        end

        for (p,w) in walk_motifs
            if !haskey(motifs, p)
                motifs[p] = 0.0
                motif_edges[p] = Set()
            end

            motifs[p] += weights[widx]
            push!(motif_edges[p], (w,widx))
        end
        widx += 1
    end
    return motifs, motif_edges
end
