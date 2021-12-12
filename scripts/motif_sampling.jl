using ArgParse, Core
include("../src/debruijnnets.jl");
include("../src/sampling.jl");
include("../src/motifs.jl")
include("../src/randomwalks.jl")
include("../src/hypa.jl")

"""
Run a hypa simulation
"""
function hypa_sample(fo, ko, ko_map, fo_map, k, num_samples, empirical_motifs)
    ko_rev_map = Dict(val=>key for (key,val) in ko_map)
    # Get the adjacency and Xi matrices
    adjacency = adjacency_matrix(ko)
    pvals, Xi = hypa(k, ko, ko_map, fo, fo_map)
    sampled_counts = Dict(m=>Dict("frequency"=>empirical_motifs[m], "samples"=>zeros(num_samples)) for m in keys(empirical_motifs))
    Threads.@threads for run in range(1, num_samples)
        # Draw a sample
        Core.println("Sampling adjacency.")
        sampled_adj = draw_sample(adjacency, Xi)
        Core.println("Sampled.")
        # Transform into walks
        cidxs = findall(>(0), sampled_adj)
        walks = Vector{Tuple}(undef, length(cidxs))
        weights = Vector{Int64}(undef, length(cidxs))
        widx = 1
        for cidx in cidxs
            u, v = Tuple(cidx)
            walk = (ko_rev_map[u]...,ko_rev_map[v][end])
            walks[widx] = walk
            weight = sampled_adj[cidx]
            weights[widx] = weight
            widx += 1
        end
        # Count motifs
        motifs = count_motifs(walks, weights)
        for (motif, val) in motifs
            sampled_counts[motif]["samples"][run] = val
        end
    end
    return sampled_counts
end

"""
Run a Biased DeBruijn graph simulation.
"""
function bdg(fo, ko, ko_map, fo_map, k, M, num_samples, empirical_motifs)
    # Get all of the walks from the korder graph
    all_walks, weights = walks_from_edges(ko, ko_map, fo_map)
    sampled_counts = Dict(m=>Dict("frequency"=>empirical_motifs[m], "samples"=>zeros(num_samples)) for m in keys(empirical_motifs))
    Threads.@threads for run in range(1, num_samples)
        # Sample M walks at a time
        walks = sample(all_walks, M)
        # Count motifs
        motifs = count_motifs(walks)
        for (motif, val) in motifs
            sampled_counts[motif]["samples"][run] = val
        end
    end
    return sampled_counts
end

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
    Run a random-walk sampling simulation.
"""
function rw(fo, k::Integer, M::Integer , num_samples::Integer, empirical_motifs,
        bias_nodes::Bool=false, weighted::Bool=false)
    nodes, A, path_counts = filter_nodes(fo, k)
    node_probabilities = compute_node_probs(nodes, path_counts)
    # Sample a bunch of times
    sampled_counts = Dict(m=>Dict("frequency"=>empirical_motifs[m], "samples"=>zeros(num_samples)) for m in keys(empirical_motifs))
    Threads.@threads for run in range(1, num_samples)
        walks = random_walks(fo, k, M, nodes, node_probabilities, bias_nodes, weighted)
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
weighted = false
if ensemble == "rw-w"
    weighted = true
end
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
elseif ensemble == "rw-uw" || ensemble == "rw-w"
    sampled_counts = rw(fo, k, M, num_samples, empirical_motifs, false, weighted)
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_e-$(ensemble).csv"
elseif ensemble == "bdg"
    sampled_counts = bdg(fo, ko, ko_map, fo_map, k, M, num_samples, empirical_motifs)
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_e-$(ensemble).csv"
elseif ensemble == "hypa"
    sampled_counts = hypa_sample(fo, ko, ko_map, fo_map, k, num_samples, empirical_motifs)
    output_filename = "../../debruijn-nets/results/motifs/$(ngram_filename)_k-$(k)_e-h.csv"
end

println(sampled_counts)
write_counts(sampled_counts, output_filename)
