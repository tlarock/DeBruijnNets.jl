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
            if length(w) != length(first(walk_edges))
                println(w)
            end
            unobserved += 1
        end
    end
    return observed, unobserved
end


"""

"""
function sample_properties(input_file, k, frequency, walk_interval, num_intervals, num_samples, output_file)
    print_interval = 100
    println("Reading graph.")
    fo, fo_map, ko, ko_map = from_ngram(input_file, frequency, k)
    rev_fo_map = Dict{Integer, String}(val=>key for (key,val) in fo_map)
    println("Done.")
    println("Converting observed kth-order edges to walks.")
    walk_edges, weights = walks_from_edges(ko, ko_map, fo_map)
    println("Done.")
    walk_edges_set = Set{Tuple}(walk_edges)
    println("Computing all possible k-edge walks.")
    all_kedge_walks, nodes = get_all_walks_with_nodes(fo, k)
    total_korder_nodes = length(nodes)
    println("Done.")
    observed_sampled = zeros(Float64, num_samples, num_intervals)
    unobserved_sampled = zeros(Float64, num_samples, num_intervals)
    missing_nodes = zeros(Float64, num_samples, num_intervals)
    missing_edges = zeros(Float64, num_samples, num_intervals)
    Threads.@threads for run in range(1, num_samples)
        println("Run: $run")
        sampled_walks = sample(all_kedge_walks, walk_interval)
        sampled_set = Set(sampled_walks)
        cmp_set = deepcopy(sampled_set)
        interval_count = 0
        intervals_run = 0
        num_observed = 0
        num_unobserved = 0
        for i in range(1, num_intervals)
            interval_count += 1
            if interval_count == print_interval
                println("i: $i, num walks: $(i*walk_interval)")
                interval_count = 0
            end

            # Construct a kth-order graph from these walks
            fo_s, ko_s, ko_map = from_walks(sampled_walks, k, rev_fo_map, ko_map)
            # Compare sampled walks to observed walks
            no, nu = compare_nodes(cmp_set, walk_edges_set)
            num_observed += no
            num_unobserved += nu
            # Update ko stats
            observed_sampled[run,i] = num_observed / ne(ko)
            unobserved_sampled[run,i] = num_unobserved / total_korder_nodes

            # Update fo stats
            missing_nodes[run,i] = nv(fo_s) / nv(fo)
            missing_edges[run,i] = ne(fo_s) / ne(fo)

            # Compute next set of walks
            nxt_walks = sample(all_kedge_walks, walk_interval)
            append!(sampled_walks, nxt_walks)
            nxt_walks_set = Set(nxt_walks)
            # Only analyze new walks
            cmp_set = setdiff(nxt_walks_set,sampled_set)
            union!(sampled_set, nxt_walks_set)
        end
    end

    x = [walk_interval*i for i in range(1, num_intervals)]
    observed_sampled_stats = mean_and_std(observed_sampled, 1)
    unobserved_sampled_stats = mean_and_std(unobserved_sampled, 1)
    missing_nodes_stats = mean_and_std(missing_nodes, 1)
    missing_edges_stats = mean_and_std(missing_edges, 1)

    open(output_file, write=true) do file
        out_str = "$(join(x, ','))\n" 
        write(file, out_str)
        for arr in [observed_sampled_stats, unobserved_sampled_stats, missing_nodes_stats, missing_edges_stats]
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
    output_file = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_samplestats.csv"
    sample_properties(input_filepath, k, frequency, walk_interval, num_intervals, runs, output_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
