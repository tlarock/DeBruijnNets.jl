mkdir -p ../results/motifs/

threads=$1

for ensemble in "cdg" "odg" "hypa" "rw-uw" "rw-w"
do
	for k in 2 3
	do
		julia --project=.. --threads $threads motif_sampling.jl -i ../data/flights/2020/coupons_2020_q1.ngram -e $ensemble -k $k -f -w -1 -o ../results/motifs/ &
		julia --project=.. --threads $threads motif_sampling.jl -i ../data/flights/2020/coupons_2020_q2.ngram -e $ensemble -k $k -f -w -1 -o ../results/motifs/

		julia --project=.. --threads $threads motif_sampling.jl -i ../data/wikispeedia/wikispeedia_finished_all.ngram -e $ensemble -k $k -w -1 -o ../results/motifs/ &
		julia --project=.. --threads $threads motif_sampling.jl -i ../data/wikispeedia/wikispeedia_unfinished_all.ngram -e $ensemble -k $k -w -1 -o ../results/motifs/
	done
done
