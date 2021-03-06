using StatsBase, Graphs, SimpleWeightedGraphs
include("sampling.jl")

"""
"""
function random_walks(G, k::Integer, num_walks::Integer, nodes, node_probabilities,
        bias_nodes::Bool=false, weighted::Bool=false)
    walks = Vector{Tuple}(undef, num_walks)
    Threads.@threads for i in range(1, num_walks)
        # Initializ walk to garbage value (always enters loop)
        walk = [0]
	while length(walk) < k+1
            start_node = sample_start_node(nodes, node_probabilities, bias_nodes)
            # Start node must have positive out degree
            while sum(G.weights[start_node,:]) == 0
                start_node = sample_start_node(nodes, node_probabilities, bias_nodes)
	    end
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
function edge_probabilities(G, node, nbrs)
    edge_probs = convert(Vector{Float64}, G.weights[node,nbrs])
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
    push!(walk, start_node)
    curr_node = start_node
    max_attempts = 5
    attempts = 0
    while length(walk) < len
        nbrs = findall(>(0), G.weights[curr_node,:])
        edge_probs = edge_probabilities(G, curr_node, nbrs)
        nxt_node = sample(nbrs, edge_probs, 1)[1]
        tries = 0
        mtries = length(nbrs)*2
        while sum(G.weights[nxt_node,:]) == 0 && tries < mtries
            nxt_node = sample(nbrs, edge_probs, 1)[1]
            tries += 1
        end

        if tries < mtries
            push!(walk, nxt_node)
            curr_node = nxt_node
        else
            # start over if we get stuck
            walk = Vector{Int64}()
            push!(walk, start_node)
            curr_node = start_node
            attempts += 1
            if attempts > max_attempts
                break
            end
        end
    end
    return Tuple(walk)
end
