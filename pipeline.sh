#!/bin/bash

mkdir -p data/split

RM="myzusPersicae_NS-families.fa"
GENOME="NS.draft.final.reviewed.sealed.fasta"
SPECIES="Myzus_persicae_NS"
export GENOME

# Cluster RepeatModeler output
cd-hit-est -n 10 -c 0.95 -i data/${RM} -o data/clustered_${RM}

# Splt RepeatModeler library into single sequences and create list of said sequences
Rscript splitter.R -t nt -f data/clustered_${RM} -o data/split/ -p 128
ls data/split/ | sed 's/.*\///' > data/queries.txt

# Search genome for each sequence and compile results
parallel --env GENOME --bar --jobs 8 -a data/queries.txt blastn -task dc-megablast -query data/split/{} -db seq/$GENOME -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out out/initial_search/{}.out -num_threads 1
cat out/initial_search/clustered_${RM}_seq_*.fasta.out > out/${GENOME}_initial_search.out

# Get multifasta from blast output ready for self alignments (adjust flanks if necessary in re-runs)
Rscript self_blast_setup.R  -g $GENOME -s $SPECIES -r $RM
ls out/initial_seq/${GENOME}*.fasta | sed 's/.*\///' > out/${GENOME}_self_queries.txt

# Local align multifastas to self ready for alignment
parallel --env GENOME --bar --jobs 8 -a out/$GENOME"_self_queries.txt" blastn -task dc-megablast -query out/initial_seq/{} -subject out/initial_seq/{} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out out/self_search/{}.out -num_threads 1

# Perform preliminary trim based on self blast
Rscript mafft_setup.R -g $GENOME -s $SPECIES

# align sequences
while read a; do
  mafft --thread 8 --localpair --adjustdirectionaccurately out/to_align/$a > out/mafft/$a;
done<"out/"$GENOME"_self_queries.txt"
