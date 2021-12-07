using StatsBase, Graphs, SimpleWeightedGraphs
include("sampling.jl")

"""
"""
function random_walks(G, k::Integer, num_walks::Integer, nodes, node_probabilities, bias_nodes::Bool=false)
    walks = Vector{Tuple}(undef, num_walks)
    Threads.@threads for i in range(1, num_walks)
        start_node = sample_start_node(nodes, node_probabilities, bias_nodes)
        walk = Tuple(randomwalk(G, start_node, k+1))
        while length(walk) < k+1
            walk = Tuple(randomwalk(G, start_node, k+1))
        end
        walks[i] = walk
    end
    return walks
end

