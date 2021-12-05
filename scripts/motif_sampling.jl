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

num_samples = 10
num_cpus = Threads.nthreads()
input_filepath = "../../debruijn-nets/data/flights/old_data/coupons_2000_01_tiny.ngram"
splitpath = split(input_filepath, '/')
ngram_filename = splitpath[end]

for k in [2]
    println("k: $k")
    fo, fo_map, ko, ko_map = from_ngram(input_filepath, true, k);
    # Count motifs
    walks, walk_freqs = walks_from_edges(ko, ko_map, fo_map)
    empirical_motifs = count_motifs(walks, walk_freqs)
    println(empirical_motifs)
    M = sum(values(empirical_motifs))
    Threads.@threads for walks_per_node in [-1, 1, 10, 50, 100]
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
	    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_ensemble-lg-s_bn_ww_wpn-$(walks_per_node).csv"
	else
	    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_ensemble-lg-s_bn_ww_wpn-$(walks_per_node).csv"
	end
	write_counts(sampled_counts, output_filename)
    end
end


