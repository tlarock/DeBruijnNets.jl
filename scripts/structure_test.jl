using ArgParse, Core, StatsBase
include("../src/debruijnnets.jl");
include("../src/sampling.jl");
include("../src/motifs.jl")
include("../src/randomwalks.jl")
include("../src/hypa.jl")

function compare_nodes(walks, walk_edges)
    observed = 0
    unobserved = 0
    for w in walks
        if w in walk_edges
            observed += 1
        else
            unobserved += 1
        end
    end
    return observed, unobserved
end


"""
Run a line-graph sampling simulation.
"""
function sample_properties(input_file, k, frequency, walk_interval, num_intervals, num_samples, output_file)
    println("Reading graph.")
    fo, fo_map, ko, ko_map = from_ngram(input_file, frequency, k)
    println("Done.")
    walk_edges, weights = walks_from_edges(ko, ko_map, fo_map)
    walk_edges = Set{Tuple}(walk_edges)
    println("Computing all k-edge walks.")
    all_kedge_walks = get_all_walks(fo, k)
    println("Done.")
    ko_rev_map = Dict(val=>key for (key, val) in ko_map)
    observed_sampled = zeros(Int64, num_samples, num_intervals)
    unobserved_sampled = zeros(Int64, num_samples, num_intervals)
    missing_nodes = zeros(Int64, num_samples, num_intervals)
    missing_edges = zeros(Int64, num_samples, num_intervals)
    for run in range(1, num_samples)
        println("Run: $run")
        sampled_walks = sample(all_kedge_walks, walk_interval)
        for i in range(1, num_intervals)
            println("i: $i, num walks: $(i*walk_interval)")
            # Construct a kth-order graph from these walks
            fo_s, ko_s, new_ko_map = from_walks(sampled_walks, k, fo_map, ko_map)
            new_rev_ko_map = Dict(val=>key for (key,val) in new_ko_map)
            # Get the set of kth order nodes
            num_observed, num_unobserved = compare_nodes(sampled_walks, walk_edges)
            observed_sampled[run,i] = num_observed
            unobserved_sampled[run,i] = num_unobserved
            missing_nodes[run,i] = nv(fo)-nv(fo_s)
            missing_edges[run,i] = ne(fo)-ne(fo_s)
            append!(sampled_walks, sample(all_kedge_walks, walk_interval))
        end
    end

    x = [walk_interval*i for i in range(1, num_intervals)]
    observed_sampled = mean_and_std(observed_sampled, 1)
    unobserved_sampled = mean_and_std(unobserved_sampled, 1)
    missing_nodes = mean_and_std(missing_nodes, 1)
    missing_edges = mean_and_std(missing_edges, 1)

    open(output_file, "w") do file
        out_str = "$(join(x, ','))\n" 
        write(file, out_str)
        for arr in [observed_sampled, unobserved_sampled, missing_nodes, missing_edges]
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
    num_intervals = arguments["iterations"]
    num_cpus = Threads.nthreads()
    splitpath = split(input_filepath, '/')
    ngram_filename = splitpath[end]
    println("k: $k")
    println(runs)
    output_file = "../../debruijn-nets/results/motifs/$(ngram_filename)_samplestats.csv"
    sample_properties(input_filepath, k, frequency, walk_interval, num_intervals, runs, output_file)
end
