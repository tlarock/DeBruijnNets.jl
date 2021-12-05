using Graphs, SimpleWeightedGraphs, StatsBase


"""
    Compute all k-edge walks from start_node
"""
function all_walks(adjacency_dict, start_node, k)
    curr_walks = Vector{Tuple}()
    neighbors = adjacency_dict[start_node]
    for ne in neighbors
	    edge = Tuple((start_node, ne))
	    push!(curr_walks, Tuple(edge))
    end
    for i in range(2, k)
        new_walks = Vector{Tuple}()
        for walk in curr_walks
	    neighbors = adjacency_dict[walk[end]]
            if length(neighbors) > 0
                for ne in neighbors
                    push!(new_walks, (walk..., ne))
                end
            end
        end
        curr_walks = new_walks
    end
    return curr_walks
end


"""
    Compute walks_per_node k-edge random walks through G
    starting from start_node.
"""
function random_walks(G, start_node, k, walks_per_node)
    curr_walks = Vector{Tuple}()
    for i in range(1, walks_per_node)
        walk = Tuple(randomwalk(G, start_node, k+1))
        while length(walk) < k+1
            walk = randomwalk(G, start_node, k)
        end
        push!(curr_walks, walk)
    end
    return curr_walks
end

    
"""
    Get the out-degree weighted probability of each walk
"""
function get_walk_probabilities!(walks, walk_weight_lookup, out_degrees)
    walk_weights = Vector{Float64}()
    for w in walks
        if !haskey(walk_weight_lookup, w)
            prod = out_degrees[w[1]]
            for i in range(2, length(w)-1)
                prod *= out_degrees[i]
            end
            walk_weight_lookup[w] = prod
        end
        push!(walk_weights, walk_weight_lookup[w])
    end
    walk_weights /= sum(walk_weights)
    return walk_weights
    end
    

"""
    Sample from the complete DeBruijn graph defined by G.
"""
function uniform_walk_sample(G::SimpleWeightedDiGraph, k::Integer, num_walks::Integer,
        walks_per_node::Integer, weight_walks::Bool, bias_nodes::Bool)
    # Compute path counts from adjacency matrix
    A = adjacency_matrix(G)
    A[findall(>(0), A)] .= 1
    # dims=2 is row sum (out walks)
    out_degrees = sum(A, dims=2)
    path_counts = sum(A^k, dims=2)
    adjacency_dict = Dict()
    for i in range(1, nv(G))
        adjacency_dict[i] = findall(>(0), A[i,:])
    end
    node_weights = Vector()
    nodes = range(1, length(path_counts))
    for node in nodes
	push!(node_weights, path_counts[node])
    end
    node_probabilities = ProbabilityWeights(node_weights / sum(node_weights))
    
    walks = Vector{Tuple}()
    walks_by_node = Dict()
    walk_weight_lookup = Dict()
    while length(walks) < num_walks
        if bias_nodes
            start_node = sample(nodes, node_probabilities)[1]
        else
            start_node = rand(nodes)[1]
        end
        
        if walks_per_node > 0
            if !haskey(walks_by_node, start_node)
                curr_walks = random_walks(G, start_node, k, walks_per_node)
                walks_by_node[start_node] = curr_walks
            else
                if length(walks_by_node) < path_counts[start_node]
                    curr_walks = random_walks(G, start_node, k, walks_per_node)
                    push!(walks_by_node[start_node], curr_walks)
                else
                    curr_walks = walks_by_node[start_node]
                end
            end
        else
            curr_walks = all_walks(adjacency_dict, start_node, k)
            walks_by_node[start_node] = curr_walks
        end
        
        # Sample one walk from walks_by_node
        if weight_walks
            walk_probabilities = get_walk_probabilities!(curr_walks, walk_weight_lookup, out_degrees)
        else
            walk_probabilities = Vector{Float64}(range(1, length(curr_walks)))
            fill!(walk_probabilities, 1. / length(curr_walks))
        end
       	if length(walk_probabilities) > 0 
        	walk = sample(curr_walks, ProbabilityWeights(walk_probabilities))
        	push!(walks, walk)
	end
    end
    return walks
end
