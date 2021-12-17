using Graphs, SimpleWeightedGraphs, StatsBase


"""
"""
function get_adj_dict(G, A)
    adjacency_dict = Dict()
    for i in range(1, nv(G))
        adjacency_dict[i] = findall(>(0), A[i,:])
    end
    return adjacency_dict
end

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
function random_walks(G, start_node::Integer, k::Integer, walks_per_node::Integer, verbose::Bool=false)
    curr_walks = Vector{Tuple}()
    max_tries = walks_per_node*4
    for i in range(1, walks_per_node)
	tries = 0
        walk = Tuple(randomwalk(G, start_node, k+1))
        while length(walk) < k+1
            walk = randomwalk(G, start_node, k+1)
	    tries += 1
	    if tries == max_tries && verbose
		println("Couldn't find a walk in $max_tries tries. \
			Returning nothing.")
		return nothing
	    end
        end
	push!(curr_walks, Tuple(walk))
    end
    return curr_walks
end

    
"""
    Get the out-degree weighted probability of each walk
"""
function get_walk_probabilities!(walks, walk_weight_lookup, out_degrees, weight_walks)
    if weight_walks
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
    else
        walk_probabilities = Vector{Float64}(range(1, length(walks)))
        fill!(walk_probabilities, 1. / length(walks))
    end
    return walk_weights
    end

function sample_start_node(nodes, node_probabilities, bias_nodes)
    if bias_nodes
        start_node = sample(nodes, node_probabilities)[1]
    else
        start_node = rand(nodes)[1]
    end
    return start_node
end


"""
"""
function filter_nodes(G, k)
    # Compute path counts from adjacency matrix
    A = adjacency_matrix(G)
    A_bin = adjacency_matrix(G)
    A_bin[findall(>(0), A)] .= 1
    path_counts = sum(A^k, dims=2)
    # We are only interested in nodes who have k-edge walks
    nodes = [u for u in range(1, nv(G)) if path_counts[u] > 0]
    println("Total nodes: $(nv(G)) nodes with k-edge walks $(length(nodes))")
    return nodes, A, A_bin, path_counts
end


function compute_node_probs(nodes, path_counts)
    node_weights = Vector()
    for node in nodes
	push!(node_weights, path_counts[node])
    end
    node_probabilities = ProbabilityWeights(node_weights / sum(node_weights))
    return node_probabilities
end

"""
Helper function for uniform walk sampler.
Mutates walks_by_node.
"""
function update_rw_sample!(G, start_node, k, walks_per_node, walks_by_node, path_counts)
    if !haskey(walks_by_node, start_node)
        # Compute the first set of walks from start_node
	curr_walks = random_walks(G, start_node, k, walks_per_node)
	while curr_walks == nothing
	    curr_walks = random_walks(G, start_node, k, walks_per_node)
	end
    	walks_by_node[start_node] = curr_walks
    else
        # add more walks
	if length(walks_by_node[start_node]) < path_counts[start_node]
	    curr_walks = random_walks(G, start_node, k, walks_per_node)
            while curr_walks == nothing
                curr_walks = random_walks(G, start_node, k, walks_per_node)
            end
            # Only add new walks
            curr_walks = [w for w in curr_walks if !(w in walks_by_node[start_node])] 
            if length(curr_walks) > 0
	    	append!(walks_by_node[start_node], curr_walks)
            end
        end
    end
    curr_walks = walks_by_node[start_node]
    return curr_walks
end


"""
    Sample from the complete DeBruijn graph defined by G.
"""
function uniform_walk_sample(G::SimpleWeightedDiGraph, k::Integer, num_walks::Integer,
        walks_per_node::Integer, weight_walks::Bool, bias_nodes::Bool, verbose::Bool=false)
    # Filter nodes based on whether they
    # originate any k-edge walks
    nodes, A, A_bin, path_counts = filter_nodes(G, k)
    adjacency_dict = get_adj_dict(G, A)
    # dims=2 is row sum (out deg/walks)
    out_degrees = sum(A_bin, dims=2)
    # Compute weights for nodes
    node_probabilities = compute_node_probs(nodes, path_counts)

    walks = Vector{Tuple}(undef, num_walks)
    walks_by_node = Dict()
    walk_weight_lookup = Dict()
    for w in range(1, num_walks)
	start_node = sample_start_node(nodes, node_probabilities, bias_nodes)
        if walks_per_node > 0
	    curr_walks = update_rw_sample!(G, start_node, k, walks_per_node, walks_by_node, path_counts)
        else
            curr_walks = all_walks(adjacency_dict, start_node, k)
            walks_by_node[start_node] = curr_walks
        end
       
        if length(curr_walks) > 0 
            # Sample one walk from walks_by_node
	    walk_probabilities = get_walk_probabilities!(curr_walks, walk_weight_lookup,
                                                     out_degrees, weight_walks)

            walk = sample(curr_walks, ProbabilityWeights(walk_probabilities))
            walks[w] = walk
	end
    end
    if verbose
	println("Computed $(length(walks)) length-$k walks out of $num_walks requested.")
    end
    return walks
end
