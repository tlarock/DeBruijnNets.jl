using Graphs, SimpleWeightedGraphs
include("debruijnnets.jl")

"""
Sample from a directed GNP ensemble.
"""
function gnp_graph(N, p)
    edgedict = Dict{Tuple, Integer}()
    for u in range(1, N)
        for v in range(1, N)
            if u == v
                continue
            end
            if rand() < p
                edgedict[(u, v)] = 1.
            end
        end
    end
    src, dest, weights = vectors_from_edgedict(edgedict)
    G = SimpleWeightedDiGraph(src, dest, weights)
    return G
end


"""
Sample from a directed Preferential Attachment ensemble
where a back link is sampled with uniform probability.
"""
function pa_graph(N, m)
    edgedict = Dict{Tuple, Integer}()
    inneighbors = Dict(node=> 0 for node in range(1, N))
    # Start with a clique on m nodes
    for i in range(1, m)
        for j in range(1, m)
            if i !=j
                edgedict[(i,j)] = 1
                inneighbors[j] += 1
            end
        end
    end

    # Grow the network up to N nodes
    for n in range(m+1, N)
        # Compute degree weights of nodes so far
        nodes = [u for u in range(1, n-1)]
        node_weights = [Float64(inneighbors[ne]) for ne in nodes]
        node_weights ./= sum(node_weights)
        node_weights = ProbabilityWeights(node_weights)
        # Sample neighbors for node n
        out_neighbors = sample(nodes, node_weights, m)
        for ne in out_neighbors
            edgedict[(n, ne)] = 1
            inneighbors[ne] += 1
            # With uniform probability, add a back link
            if rand() < (1. / length(out_neighbors))
                edgedict[(ne,n)] = 1
                inneighbors[n] += 1
            end
        end
    end
    src, dest, weights = vectors_from_edgedict(edgedict)
    G = SimpleWeightedDiGraph(src, dest, weights)
    return G
end
