using StatsBase, Graphs, SimpleWeightedGraphs
include("sampling.jl")

"""
"""
function random_walks(G, k::Integer, num_walks::Integer, nodes, node_probabilities,
        bias_nodes::Bool=false, weighted::Bool=false)
    walks = Vector{Tuple}(undef, num_walks)
    Threads.@threads for i in range(1, num_walks)
        start_node = sample_start_node(nodes, node_probabilities, bias_nodes)
        walk = Tuple(randomwalk(G, start_node, k+1))
        while length(walk) < k+1
            if !weighted
                walk = Tuple(randomwalk(G, start_node, k+1))
            else
                walk = weighted_rw(G, start_node, k+1)
            end
        end
        walks[i] = walk
    end
    return walks
end


"""
Compute out weight basd edge probabilitis.
"""
function edge_probabilities(G, node)
    edge_probs = G.weights[node,:]
    edge_probs ./= sum(edge_probs)
    return ProbabilityWeights(edge_probs)
end


"""
Computed a weighted random walk of length
len through G from start_node.

WARNING: Assume start_node has at least 1 len-edge
walk originating from it
"""
function weighted_rw(G, start_node, len)
    walk = Vector{Int64}()
    curr_node = start_node
    while length(walk) < len
        nbrs = outneighbors(G, curr_node)
        edge_probs = edge_probabilities(G, node)
        nxt_node = sample(nbrs, edge_probs, 1)
        tries = 0
        mtries = length(nbrs)*2
        while sum(G.weights[nxt_node,:]) == 0 && tries < mtries
            nxt_node = sample(nbrs, edge_probs, 1)
            tries += 1
        end
        if tries <= mtries
            push!(walk, nxt_node)
        else
            # start over if we get stuck
            walk = Vector{Int64}()
        end
    end
    return Tuple(walk)
end
