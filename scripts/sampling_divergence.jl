using ArgParse, StatsBase
include("../src/debruijnnets.jl");
include("../src/sampling.jl");
include("../src/motifs.jl");
include("../src/randomgraphs.jl");
include("../src/randomwalks.jl");

function empty_dict(counts)
    return Dict(m=>0 for m in keys(counts))
end

function divergence_by_nodes(G, k, true_counts, true_dist, walk_interval,
                              num_intervals,  all_kwalks, methods)

    nodes, A, path_counts = filter_nodes(G, k)
    node_probabilities = compute_node_probs(nodes, path_counts)
    sampled_counts = Dict(key=>empty_dict(true_counts) for key in methods)
    divergences = Dict(key=>Vector{Float64}(undef, num_intervals) for key in keys(sampled_counts))
    for i in range(1, num_intervals)
        for met in methods
            if met == "all"
                walks = uniform_walk_sample(G, k, walk_interval,
                                            -1, true, true)
	    elseif met == "true"
		walks = sample(all_kwalks, walk_interval)
            elseif met == "rw"
                walks = random_walks(G, k, walk_interval, nodes, node_probabilities, false, false)
            else
                walks = uniform_walk_sample(G, k, walk_interval,
                                            met, true, true)
            end

            curr_counts = count_motifs(walks, ones(length(walks)))
            for (motif, val) in curr_counts
                sampled_counts[met][motif] = val
            end

            sampled_dist = [Float64(sampled_counts[met][motif]) for motif in keys(true_counts)]
            sampled_dist .+= 10^-12
            sampled_dist ./= sum(sampled_dist)
            div = kldivergence(true_dist, sampled_dist)
            divergences[met][i] = div
        end
    end
    return divergences
end

function write_to_file(output, output_file)
    open(output_file, "w") do file
        for (k, kdict) in output
	    for (m, mtup) in kdict
		out_str = "$k|$(m)|$(join(mtup[1], ","))|$(join(mtup[2], ","))\n"
                write(file, out_str)
            end
        end
    end
    return nothing
end

function get_truth(G, adj, k)
    adjacency_dict = get_adj_dict(G, adj)
    walklist = Vector{Tuple}()
    for node in range(1, nv(G))
        nodewalks = all_walks(adjacency_dict, node, k)
        append!(walklist, nodewalks)
    end
    true_counts = count_motifs(walklist)
    true_dist = [true_counts[motif] for motif in keys(true_counts)]
    true_dist ./= sum(true_dist)
    return true_dist, true_counts
end


function get_all_truths(G, kvals)
    true_dist_dict = Dict()
    true_count_dict = Dict()
    all_kedge_walks = Dict()
    total_walks = Dict()
    adj = adjacency_matrix(G)
    for k in kvals
        true_dist_dict[k], true_count_dict[k] = get_truth(G, adj, k)
        total_walks[k] = sum(adj^k)
        all_kedge_walks[k] = get_all_walks(G, k)
    end
    return true_dist_dict, true_count_dict, total_walks, all_kedge_walks
end

function graph_and_motifs(er, p, pa, m, kvals)
   if er
        G = gnp_graph(N, p)
    elseif pa
        G = pa_graph(N, m)
   end
   return G, get_all_truths(G, kvals)
end

function run_gpr(er::Bool, p::Float64, pa::Bool, m::Integer,
                      kvals, methods, num_runs::Integer)
    output = Dict(k=>Dict(method=>Array{Float64}(undef, num_runs, iterations) for method in methods) for k in kvals)
    Threads.@threads for i in range(1, num_runs)
        # Create a new graph for each run 
        G, (true_dist_dict, true_count_dict, total_walks, all_kedge_walks) = graph_and_motifs(er, p, pa, m, kvals)
        Threads.@threads for k in kvals
            divergences = divergence_by_nodes(G, k, true_count_dict[k], true_dist_dict[k], walk_interval, iterations, all_kedge_walks[k], methods)
            for method in keys(divergences)
                for j in range(1, length(divergences[method]))
                    output[k][method][i,j] = divergences[method][j]
                end
            end
        end
        println(i)
     end
   return output
end

function run_oneg(er::Bool, p::Float64, pa::Bool, m::Integer,
                      kvals, methods, num_runs::Integer)
    output = Dict(k=>Dict(method=>Array{Float64}(undef, num_runs, iterations) for method in methods) for k in kvals)
    # Create a single graph
    G, true_dist_dict, true_count_dict, total_walks, all_kedge_walks = G, true_dist_dict, true_count_dict, total_walks, all_kedge_walks
end
    println("WHAT THE FUCK IS HAPPENING")
    # Sample from it many times
    Threads.@threads for i in range(1, num_runs)
        Threads.@threads for k in kvals
           divergences = divergence_by_nodes(G, k, true_count_dict[k], true_dist_dict[k], walk_interval, iterations, all_kedge_walks[k], methods)
            for method in keys(divergences)
                for j in range(1, length(divergences[method]))
                   output[k][method][i,j] = divergences[method][j]
                end
            end
        end
        println(i)
    end
    return output
end


function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
        "-k"
	    help = "Order of the motifs."
	    arg_type = Int
        "-N", "--num_nodes"
	    help = "Number of nodes"
	    arg_type = Int
	"--pa"
	    action = :store_true
	    help = "Preferential attachment model."
	"-m"
	    arg_type = Int
	    default = 5
	    help = "Preferential attachment param.
                      Ignored for -e."
	"-e", "--er"
	    action = :store_true
	    help = "ER model."
	"-p"
	    arg_type = Float64
	    help = "ER model param. Ignored for -b."
	"-w", "--walk_interval"
	    help = "Number of walks per iteration."
	    arg_type = Int
	"-i", "--iterations"
	    help = "Number of iterations."
	    arg_type = Int
	"-r", "--runs"
	    help = "Number of runs."
	    arg_type = Int
	"--gpr"
	    help = "Compute a new graph every run."
	    action = :store_true
end
return parse_args(args, s)
end

function main()
	arguments = parse_commandline(ARGS)
	println(arguments)
	k = arguments["k"]
	N = arguments["num_nodes"]
	pa = arguments["pa"]
	m = arguments["m"]
	er = arguments["er"]
	p = arguments["p"]
	walk_interval = arguments["walk_interval"]
	iterations = arguments["iterations"]
	global runs = arguments["runs"]
	gpr = arguments["gpr"]
	num_cpus = Threads.nthreads()
	kvals = [2, 3, 4]
	methods = ["true", "all", "rw", 1, 50, 100]

	if !gpr
	    println(er, p, pa, m, kvals, methods, runs)
	    #output = run_oneg(er, p, pa, m, kvals, methods, runs)
	else
	    #output = run_gpr(er, p, pa, m, kvals, methods, runs)
            println(er, p, pa, m, kvals, methods, runs)

	end


	mean_output = Dict(k=>Dict() for k in keys(output))
	for (k, kdict) in output
	    for (method, method_arr) in kdict
		    mean_output[k][method] = mean_and_std(method_arr, 1)
	    end
	end

	if er
	    model = "er"
	    param = p
	else
	    model = "pa"
	    param = m
	end
	 
	output_filename = "../../debruijn-nets/results/motifs/$(model)_N-$(N)_param-$(param)_interval-$(walk_interval)_iters-$(iterations)_runs-$(runs)_kl"
	if gpr
	    output_filename = output_filename + "_gpr.csv"
	else
	    output_filename = output_filename + "_1g.csv"
	end

	write_to_file(mean_output, output_filename)
end

main()
