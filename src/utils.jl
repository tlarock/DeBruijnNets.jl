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

function get_all_walks(G, k)
    adjacency_dict = get_adj_dict(G, adjacency_matrix(G))
    walklist = Vector{Tuple}()
    for node in range(1, nv(G))
        nodewalks = all_walks(adjacency_dict, node, k)
        append!(walklist, nodewalks)
    end
    return walklist
end

