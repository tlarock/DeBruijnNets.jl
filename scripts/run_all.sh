mkdir -p ../results/motifs/

# Only argument to the script is the number of 
# threads for julia to use. Default is 1.
threads=1
if [ "$#" = 1 ]; then
	threads=$1
fi

for ensemble in "cdg" "odg" "hypa" "rw-uw" "rw-w"
do
	for k in 2 3
	do
		echo julia --project=.. --threads $threads motif_sampling.jl -i ../data/flights/2020/coupons_2020_q1.ngram -e $ensemble -k $k -f -w -1 -o ../results/motifs/ &
		echo julia --project=.. --threads $threads motif_sampling.jl -i ../data/flights/2020/coupons_2020_q2.ngram -e $ensemble -k $k -f -w -1 -o ../results/motifs/

		echo julia --project=.. --threads $threads motif_sampling.jl -i ../data/wikispeedia/wikispeedia_finished_all.ngram -e $ensemble -k $k -w -1 -o ../results/motifs/ &
		echo julia --project=.. --threads $threads motif_sampling.jl -i ../data/wikispeedia/wikispeedia_unfinished_all.ngram -e $ensemble -k $k -w -1 -o ../results/motifs/
	done
done
