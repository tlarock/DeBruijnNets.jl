using ArgParse, Core, StatsBase
include("../src/debruijnnets.jl");
include("../src/sampling.jl");
include("../src/motifs.jl")
include("../src/randomwalks.jl")
include("../src/hypa.jl")

function compare_nodes(observed_nodes, sampled_nodes)
    observed = length(intersect(observed_nodes, sampled_nodes))
    unobserved = length(sampled_nodes) - observed
    return observed, unobserved
end

function compare_edges(observed_edges, sampled_edges)
    observed = length(intersect(observed_edges, sampled_edges))
    unobserved = length(sampled_edges) - observed
    return observed, unobserved
end


"""

"""
function sample_properties(input_file, k, frequency, walk_interval, max_intervals, num_samples, output_file)
    print_interval = 100
    println("Reading graph.")
    fo, fo_map, ko, ko_map = from_ngram(input_file, frequency, k)
    rev_fo_map = Dict{Integer, String}(val=>key for (key,val) in fo_map)
    # Observed first-order nodes/edges
    observed_fo_nodes = length(fo_map)
    observed_fo_edges = ne(fo)
    # Observed kth-order nodes
    observed_nodes = Set{Tuple}(keys(ko_map))
    observed_korder_nodes = length(observed_nodes)
    println("Done.")
    println("Computing all possible k-edge walks.")
    all_kedge_walks, all_korder_nodes = get_all_walks_with_nodes(fo, k)
    total_korder_nodes = length(all_korder_nodes)
    total_korder_edges = length(all_kedge_walks)
    println("Done.")
    println("observed_korder_nodes: $(observed_korder_nodes)")
    println("total_korder_nodes: $(total_korder_nodes)")
    observed_edges, weights = walks_from_edges(ko, ko_map, fo_map)
    observed_edgeset = Set{Tuple}(observed_edges)
    observed_ko_edges = length(observed_edges)
    observed_sampled_nodes = zeros(Float64, num_samples, max_intervals)
    unobserved_sampled_nodes = zeros(Float64, num_samples, max_intervals)
    observed_sampled_edges = zeros(Float64, num_samples, max_intervals)
    unobserved_sampled_edges = zeros(Float64, num_samples, max_intervals)
    missing_nodes = zeros(Float64, num_samples, max_intervals)
    missing_edges = zeros(Float64, num_samples, max_intervals)
    Threads.@threads for run in range(1, num_samples)
        println("Run: $run")
        sampled_walks = sample(all_kedge_walks, walk_interval)
        sampled_edgeset = Set{Tuple}(sampled_walks)
        interval_count = 0
        intervals_run = 0
        # Redefining makes this type stable
        non::Integer = 0
        nun::Integer = 0
        noe::Integer = 0
        nue::Integer = 0
        i = 1
        while (nun / total_korder_nodes) < 1.0 && i <= max_intervals
            # Construct a kth-order graph from these walks
            fo_s, ko_s, new_ko_map = from_walks(sampled_walks, k, rev_fo_map)
            sampled_nodes = Set{Tuple}(keys(new_ko_map))
            
            # Compare nodes/edges from sampled walks to nodes from observed walks
            non, nun = compare_nodes(observed_nodes, sampled_nodes)
            noe, nue = compare_edges(observed_edgeset, sampled_edgeset)

            # Update ko stats
            observed_sampled_nodes[run,i] = non / observed_korder_nodes
            unobserved_sampled_nodes[run,i] = nun / total_korder_nodes
            observed_sampled_edges[run,i] = noe / observed_ko_edges
            unobserved_sampled_edges[run,i] = nue / total_korder_edges

            # Update fo stats
            missing_nodes[run,i] = nv(fo_s) / observed_fo_nodes
            missing_edges[run,i] = ne(fo_s) / observed_fo_edges

            # Compute next set of walks
            nxt_walks = sample(all_kedge_walks, walk_interval)
            append!(sampled_walks, nxt_walks)
            union!(sampled_edgeset, Set(nxt_walks))
            # Print for logging
            interval_count += 1
            if interval_count == print_interval
                println("i: $i, num walks: $(i*walk_interval), proportion unobserved: $(nun/total_korder_nodes)")
                interval_count = 0
            end
            i += 1
        end
    end

    x = [walk_interval*i for i in range(1, max_intervals)]
    observed_sampled_nodes_stats = mean_and_std(observed_sampled_nodes, 1)
    unobserved_sampled_nodes_stats = mean_and_std(unobserved_sampled_nodes, 1)
    observed_sampled_edges_stats = mean_and_std(observed_sampled_edges, 1)
    unobserved_sampled_edges_stats = mean_and_std(unobserved_sampled_edges, 1)

    missing_nodes_stats = mean_and_std(missing_nodes, 1)
    missing_edges_stats = mean_and_std(missing_edges, 1)

    open(output_file, write=true) do file
        out_str = "$(join(x, ','))\n" 
        write(file, out_str)
        for arr in [observed_sampled_nodes_stats, unobserved_sampled_nodes_stats, observed_sampled_edges_stats, unobserved_sampled_edges_stats, missing_nodes_stats, missing_edges_stats]
            out_str = "$(join(arr[1], ','))|$(join(arr[2], ','))\n" 
            write(file, out_str)
        end
    end
end

function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
	"--file", "-g"
            help = "Full path to input file"
	"--frequency", "-f"
	    help = "Frequencies in ngram?"
	    action = :store_true
        "-k"
	    help = "Order of the motifs."
	    arg_type = Int
    	"-w", "--walk_interval"
	    help = "Number of walks per iteration."
	    arg_type = Int
	"-i", "--iterations"
	    help = "Number of iterations."
	    arg_type = Int
	"-r", "--runs"
	    help = "Number of runs."
	    arg_type = Int
    end
return parse_args(args, s)
end

function main()
    arguments = parse_commandline(ARGS)
    println(arguments)
    input_filepath = arguments["file"]
    println(input_filepath)
    frequency = arguments["frequency"]
    println(frequency)
    k = arguments["k"]
    runs = arguments["runs"]
    walk_interval = arguments["walk_interval"]
    max_intervals = arguments["iterations"]
    num_cpus = Threads.nthreads()
    splitpath = split(input_filepath, '/')
    ngram_filename = splitpath[end]
    println("k: $k")
    println(runs)
    output_file = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_samplestats.csv"
    sample_properties(input_filepath, k, frequency, walk_interval, max_intervals, runs, output_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
