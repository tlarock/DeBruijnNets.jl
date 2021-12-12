using Graphs, SimpleWeightedGraphs, StatsBase, SparseArrays, LinearAlgebra, Distributions


"""
Given a DeBruijn graph ko, its corresponding
first order graph fo, and an integer order k,
construct the Xi Matrix for ko and compute
p-values on every edge.
"""
function hypa(k::Integer,
        ko::SimpleWeightedDiGraph, ko_map::Dict{Tuple, Integer},
        fo::SimpleWeightedDiGraph, fo_map::Dict{String, Integer})
    # Construct Xi
    Xi = compute_xi(ko, ko_map, fo, fo_map)
    pvals = compute_pvals(adjacency_matrix(ko), Xi, false) 
    return pvals
end

"""
Compute p-values
"""
function compute_pvals(adjacency, Xi, log_p::Bool=false)
    total_observations = sum(adjacency)
    total_xi = sum(Xi)
    nonzero_xi = findall(>(0), Xi)
    pvals = spzeros(Xi.m, Xi.n)
    for cidx in nonzero_xi
        observed = adjacency[cidx]
        xi = Xi[cidx]
        hy  = Hypergeometric(total_observations, total_xi-total_observations, xi)
        if log_p
            pvals[cidx] = logcdf(hy, observed)
        else
            pvals[cidx] = cdf(hy, observed)
        end
    end
    return pvals
end

"""
Draw a sample from the ensemble defined
by adjacency and Xi.
"""
function draw_sample(adjacency, Xi)
    sampled_adjacency = spzeros(adjacency.m, adjacency.n)
    xi_sum = sum(Xi)
    adj_sum = sum(adjacency)
    xi_accum = 0
    for cidx in findall(>(0), Xi)
        xi_accum += Xi[cidx]
        hy = Hypergeometric(Xi[cidx], xi_sum-xi_accum, adj_sum)
        sampled_adjacency[cidx] = rand(hy)
        adj_sum -= sampled_adjacency[cidx]
    end
    return sampled_adjacency
end

"""
FUNCTIONS FOR COMPUTING THE XI MATRIX
"""

"""
Compute Xi matrix.
"""
function compute_xi(
        ko::SimpleWeightedDiGraph, ko_map::Dict{Tuple, Integer},
        fo::SimpleWeightedDiGraph, fo_map::Dict{String, Integer},
        tol::Float64=0.0001)
    fo_adj = adjacency_matrix(fo)
    ko_adj = adjacency_matrix(ko)
    fo_rev_map = Dict(val=>key for (key,val) in fo_map)
    ko_rev_map = Dict(val=>key for (key,val) in ko_map)
    Xi = spzeros(Float64, nv(ko), nv(ko))
    for ko_node in range(1, nv(ko))
        # Get last first-order node in node
        ko_node_mapped = ko_rev_map[ko_node]
        last_node = ko_node_mapped[end]
        fo_node = fo_map[last_node]
        # Get all neighbors of last node
        neighbors = findall(>(0), fo_adj[fo_node,:])
        # for each neighbor, check if corresponding
        # higher-order node is in ko
        prefix_mapped = [u for u in ko_node_mapped[2:end]]
        out_deg = sum(ko_adj[ko_node, findall(>(0), ko_adj[ko_node,:])])
        if out_deg < 1
            continue
        end
        for ne in neighbors
            candidate_node = Vector(prefix_mapped)
            push!(candidate_node, fo_rev_map[ne])
            candidate_node = Tuple([u for u in candidate_node])
            if haskey(ko_map, candidate_node) 
                # if so, compute a xival. Else do nothing.
                candidate_idx = ko_map[candidate_node]
                in_deg = sum(ko_adj[findall(>(0), ko_adj[:,candidate_idx]), candidate_idx])
                Xi[ko_node, candidate_idx] = out_deg * in_deg
            end
        end
    end
    fit_xi!(Xi, ko, tol)
    return Xi
end

"""
Fit the Xi matrix by redistributing
the excess degree from impossible
edges.
"""
function fit_xi!(Xi, ko, tol::Float64)
    adj = adjacency_matrix(ko)
    m = sum(adj)
    indegs = sum(adj, dims=1)
    outdegs = sum(adj, dims=2)
    xi_sum = sum(Xi)
    rmse = compute_rmse(indegs, outdegs, Xi, xi_sum, m)
    println("rmse: $rmse") 
    while rmse > tol
        # Fit rows
        Xi = xifit_row(m, Xi, outdegs)
        # Fit columns
        Xi = transpose(xifit_row(m, transpose(Xi), indegs))

        xi_sum = sum(Xi)
        rmse_new = compute_rmse(indegs, outdegs, Xi, xi_sum, m)

        if rmse <= rmse_new || rmse < tol
            break
        end
        rmse = rmse_new
    end
    return Xi
end

"""
Operations to redistribute for a single row
or column.
"""
function xifit_row(m, Xi, degs)
    xi_sum = sum(Xi)
    Xi .= (Xi .* m^2) ./ xi_sum
    xi_sum = sum(Xi)

    exp_degs = vec(sum((Xi.*m)./xi_sum, dims=2))
    nonzero = findall(!=(0), exp_degs)
    ratio = spzeros(length(exp_degs))
    ratio[nonzero] .= degs[nonzero]./exp_degs[nonzero]
    Xi = ratio.*Xi
    return Xi
end


"""
Compute RMSE
"""
function compute_rmse(indegs, outdegs, Xi, xi_sum, m)
    in_sum = sum((Xi.*m) ./ xi_sum, dims=1)
    in_rmse = (in_sum-indegs).^2
    in_rmse = sum(sqrt.(in_rmse)./2)

    out_sum = sum((Xi.*m) ./ xi_sum, dims=2)
    out_rmse = (out_sum-outdegs).^2
    out_rmse = sum(sqrt.(out_rmse)./2)

    return out_rmse + in_rmse
end
