using ArgParse
include("../src/debruijnnets.jl");
include("../src/sampling.jl");
include("../src/motifs.jl")

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
num_cpus = Threads.nthreads()
splitpath = split(input_filepath, '/')
ngram_filename = splitpath[end]

println("k: $k")
fo, fo_map, ko, ko_map = from_ngram(input_filepath, frequency, k);
# Count motifs
walks, walk_freqs = walks_from_edges(ko, ko_map, fo_map)
empirical_motifs = count_motifs(walks, walk_freqs)
println(empirical_motifs)
M = sum(values(empirical_motifs))
println("wpn: $walks_per_node.")
sampled_counts = Dict(m=>Dict("frequency"=>empirical_motifs[m], "samples"=>zeros(num_samples)) for m in keys(empirical_motifs))
# Sample a bunch of times
Threads.@threads for run in range(1, num_samples)
    walks = uniform_walk_sample(fo, k, M, walks_per_node, true, true)
    # Count motifs
    motifs = count_motifs(walks)
    for (motif, val) in motifs
	sampled_counts[motif]["samples"][run] = val
    end
end
println(sampled_counts)
if walks_per_node > 0
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_e-lg-s_bn_ww_wpn-$(walks_per_node).csv"
else
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_ensemble-lg-s_bn_ww_wpn-$(walks_per_node).csv"
end
write_counts(sampled_counts, output_filename)

