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
    return projected
end

"""
Accept a set of walks, assumed to be of the same length.
Project each walk into a motif. If no weights are provided,
each has weight 1.
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
