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
