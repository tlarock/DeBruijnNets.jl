# DeBruijnNets.jl

Julia code for computing and sampling from DeBruijn graphs constructed from ngram data. Implements the techniques in the following two papers:


1. LaRock, T., Scholtes, I., Eliassi-Rad, T. Sequential Motifs in Observed Walks. arXiv:2112.05642 [physics.soc-ph], December 2021. [Link.](https://arxiv.org/abs/2112.05642)
2. LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T., Schweitzer, F. HYPA: Efficient Detection of Path Anomalies in Time Series Data on Networks. In Proc. of the 2020 SIAM Int’l Conf. on Data Mining, SDM’20, 2020: 460-468. [Link.](https://epubs.siam.org/doi/abs/10.1137/1.9781611976236.52)


This package does not define a DeBruijnNet type, but instead wraps `SimpleWeightedDiGraph` from [`SimpleWeightedGraphs.jl`](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl), which is in turn built on [`Graphs.jl`](https://github.com/JuliaGraphs/Graphs.jl).

Construct a DeBruijn graph as follows:
```
first_order, fo_map, kth_order, ko_map = from_ngram(input_filepath, frequency, k)
```

`input_filepath` should point to an ngram file where each line corresponds to a comma-separated ngram. If the last entry in each row is a frequency or weight for that ngram, set parameter `frequency = true`. If the rows contain just the ngrams, set `frequency = false`. `k` is an integer corresponding to the order of the DeBruijn graph that will be constructed. The `_map` variables map each node to an index; `fo_map` maps from original node (content of ngram file) to its integer index in `first_order`, and `ko_map` maps from tuple representing a kth-order node to its integer index in the graph `kth_order`. 

Prelminary implementation provided with no guarantees and limited documentation. Please report bugs as issues.
