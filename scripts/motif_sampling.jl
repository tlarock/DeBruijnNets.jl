using ArgParse
include("../src/debruijnnets.jl");
include("../src/sampling.jl");
include("../src/motifs.jl")

"""
Run a line-graph sampling simulation.
"""
function lgs(fo, k, M, num_samples,  walks_per_node, empirical_motifs)
    println("wpn: $walks_per_node.")
    if walks_per_node < 0
        println("Computing all walks for every node.")
        all_kedge_walks = get_all_walks(fo, k)
        println("Done with computation.")
    end

    # Sample a bunch of times
    sampled_counts = Dict(m=>Dict("frequency"=>empirical_motifs[m], "samples"=>zeros(num_samples)) for m in keys(empirical_motifs))
    Threads.@threads for run in range(1, num_samples)
        if walks_per_node > 0
            walks = uniform_walk_sample(fo, k, M, walks_per_node, true, true)
        else
            walks = sample(all_kedge_walks, M)
        end
        # Count motifs
        motifs = count_motifs(walks)
        for (motif, val) in motifs
            sampled_counts[motif]["samples"][run] = val
        end
    end
    return sampled_counts
end


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

"""
    Run a random-walk sampling simulation.
"""
function rw(fo, k::Integer, M::Integer , num_samples::Integer, empirical_motifs, bias_nodes::Bool=false)
    nodes, A, path_counts = filter_nodes(fo, k)
    node_probabilities = compute_node_probs(nodes, path_counts)
    # Sample a bunch of times
    sampled_counts = Dict(m=>Dict("frequency"=>empirical_motifs[m], "samples"=>zeros(num_samples)) for m in keys(empirical_motifs))
    Threads.@threads for run in range(1, num_samples)
        walks = random_walks(fo, k, M, nodes, node_probabilities, bias_nodes)
        # Count motifs
        motifs = count_motifs(walks)
        for (motif, val) in motifs
            sampled_counts[motif]["samples"][run] = val
        end
    end
    return sampled_counts
end


function write_counts(sampled_counts, output_file)
    open(output_file, "w") do file
        for (m, mdict) in sampled_counts
            out_str = join(m, ',') * '|' * string(mdict["frequency"]) * '|' * join(mdict["samples"], ",") * "\n"
            write(file, out_str)
        end
    end
    return nothing
end

function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
	"--input_file", "-i"
        help = "Full path to input file"
        #arg_type = String
	"--frequency", "-f"
	    help = "Frequencies in ngram?"
	    action = :store_true
        "-k"
	    help = "Order of the motifs."
	    arg_type = Int
	"-w", "--walks_per_node"
	    help = "Walks per node. -1 for all."
            arg_type = Int
	"-n", "--num_samples"
	    help = "Number of samples to draw."
	    arg_type = Int
	    default = 10
	"-e", "--ensemble"
	    help = "Ensemble to use."
	    arg_type = String
    end
return parse_args(args, s)
end

arguments = parse_commandline(ARGS)

println(arguments)
input_filepath = arguments["input_file"]
println(input_filepath)
frequency = arguments["frequency"]
println(frequency)
k = arguments["k"]
walks_per_node = arguments["walks_per_node"]
num_samples = arguments["num_samples"]
ensemble = arguments["ensemble"]
num_cpus = Threads.nthreads()
splitpath = split(input_filepath, '/')
ngram_filename = splitpath[end]

println("k: $k")
fo, fo_map, ko, ko_map = from_ngram(input_filepath, frequency, k);
if has_self_loops(fo)
   println("fo has selfloops..")
end
# Count motifs
walks, walk_freqs = walks_from_edges(ko, ko_map, fo_map)
empirical_motifs = count_motifs(walks, walk_freqs)
println(empirical_motifs)
M = sum(values(empirical_motifs))
if ensemble == "lg-s"
    sampled_counts = lgs(fo, k, M, num_samples, walks_per_node, empirical_motifs)
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_e-$(ensemble)_bn_ww_wpn-$(walks_per_node).csv"
elseif ensemble == "rw-uw"
    sampled_counts = rw(fo, k, M, num_samples, empirical_motifs)
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_e-$(ensemble)_bn_ww.csv"

end

println(sampled_counts)
write_counts(sampled_counts, output_filename)
