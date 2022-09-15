using CSV, Tables
include("../src/debruijnnets.jl")
include("../src/motifs.jl")
include("../src/motif_embeddings.jl")

function compute_embedding(filename, frequency, k_vals)
    empirical_motifs = Dict()
    motif_edges = Dict()
    all_walks = Vector{Tuple}()
    all_walk_freqs = Vector{Int64}()
    fo_maps = Dict()
    for k in k_vals
        println(k)
        forder, fo_map, korder, ko_map = from_ngram(filename, frequency, k)
        fo_maps[k] = fo_map
        walks, walk_freqs = walks_from_edges(korder, ko_map, fo_map)
        all_walks = [all_walks; walks]
        all_walk_freqs = [all_walk_freqs; walk_freqs]
    end

    empirical_motifs, motif_edges = motifs_with_mapping(all_walks; weights=all_walk_freqs)

    N = length(fo_maps[2])
    weighted = true
    embedding, embedding_label = positional_embedding(motif_edges, all_walks, all_walk_freqs, weighted, N)

    CSV.write(filename * "_mapping.csv", fo_maps[k_vals[1]]; header=false)
    CSV.write(filename * "_embedding.csv", Tables.table(embedding); header=embedding_label)
    return nothing
end

function compute_exclusive_embedding(filename, frequency, k_vals)
    all_walks, fo_map = walks_from_ngram(filename, frequency)
    walks = Vector{Tuple}()
    weights = Vector{Int64}()
    for (walk, weight) in all_walks
        push!(walks, walk)
        push!(weights, weight)
    end
    empirical_motifs, motif_edges = exclusive_motifs(walks, k_vals; weights=weights)
    N = length(fo_map)
    weighted = true
    embedding, embedding_label = positional_embedding(motif_edges, walks, weights, weighted, N)

    CSV.write(filename * "_exclusive_mapping.csv", fo_map; header=false)
    CSV.write(filename * "_exclusive_embedding.csv", Tables.table(embedding); header=embedding_label)
    return nothing
end


#filename = "../data/flights/2020/coupons_2020_q1.ngram"
#frequency = true
filename = "../data/wikispeedia/wikispeedia_finished_all.ngram"
frequency = false
k_vals = [2,3]
compute_exclusive_embedding(filename, frequency, k_vals)
