#!/bin/bash

mkdir -p data/split

RM="myzusPersicae_bassLabG006-families.fa"
GENOME="myzusPersicae_bassLabG006.fasta"
export GENOME

cd-hit-est -n 10 -c 0.95 -i data/${RM} -o data/clustered_${RM}

Rscript splitter.R -t nt -f data/clustered_${RM} -o data/split/ -p 128

ls data/split/ | sed 's/.*\///' > data/queries.txt

parallel --env GENOME --bar --jobs 8 -a data/queries.txt blastn -task dc-megablast -query data/split/{} -db seq/$GENOME -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out out/initial_search/{}.out -num_threads 1

cat out/initial_search/clustered_${RM}_seq_*.fasta.out > out/${GENOME}_initial_search.out

Rscript self_blast_setup.R

ls out/initial_seq/${GENOME}*.fasta | sed 's/.*\///' > out/${GENOME}_self_queries.txt

parallel --env GENOME --bar --jobs 8 -a out/$GENOME"_self_queries.txt" blastn -task dc-megablast -query out/initial_seq/{} -subject out/initial_seq/{} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out out/self_search/{}.out -num_threads 1

Rscript mafft_setup.R

while read a; do
  mafft --thread 8 --localpair --adjustdirectionaccurately out/to_align/$a > out/mafft/$a;
done<"out/"$GENOME"_self_queries.txt"
