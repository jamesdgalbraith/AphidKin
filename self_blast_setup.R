#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-g", "--genome_name"), type="character", default=NULL, 
              help="genome name", metavar="character"),
  make_option(c("-s", "--species_name"), type="character", default=NULL,
              help="species name", metavar="character"),
  make_option(c("-r", "--repeatmodeler"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-5", "--flank5"), type="character", default=NULL,
              help="Length to extend 5' flank", metavar="integer"),
  make_option(c("-3", "--flank3"), type="character", default=NULL,
              help="Length to extend 3' flank", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$genome_name)) {
  stop("Genome name is needed")
} else {
  # set genome names
  genome_name <- opt$genome_name
}

if (is.null(opt$species_name)) {
  species_name <- opt$genome_name
} else {
  species_name <- opt$species_name
}

# Set flanking region length if not manually set
if (is.null(opt$flank5)) {
  flank5 <- 2500
} else {
  flank5 <- opt$flank5
}
if (is.null(opt$flank3)) {
  flank3 <- 2500
} else {
  flank3 <- opt$flank3
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

# read in blast file
blast_out <- read_tsv(file = paste0("out/", opt$genome_name, "_initial_search.out"),
                      col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend",
                                    "qlen", "sstart", "send", "slen", "evalue", "bitscore"))

# filter blast results
blast_out_trimmed <- blast_out %>%
  dplyr::filter(
    pident >= 80,
    length >= 0.5 * qlen
    ) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:20) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-"),
         start = base::ifelse(sstart < send, sstart - flank5, send - flank3),
         end = base::ifelse(sstart > send, sstart + flank3, send + flank5),
         start = ifelse(start < 1, 1, start),
         end = ifelse(end > slen, slen, end),
         seqnames = paste0(qseqid, "___", seqnames)) %>%
  dplyr::select(seqnames, qseqid, start, end, strand) %>%
  plyranges::as_granges()

# collapse overlapping intervals
blast_out_merged <- as_tibble(plyranges::reduce_ranges(blast_out_trimmed)) %>%
  mutate(seqnames = as.character(seqnames)) %>%
  separate(seqnames, into = c("qseqid", "seqnames"), sep = "___") %>%
  dplyr::select(-strand, -width) %>%
  plyranges::as_granges()

# ensure more than two sequences, identify those with less
blast_out_tbl <- tibble(qseqid = names(table(blast_out_merged$qseqid)), n = as.integer(table(blast_out_merged$qseqid))) %>%
  filter(n > 2)
blast_out_ignored_tbl <- tibble(qseqid = names(table(blast_out_merged$qseqid)), n = as.integer(table(blast_out_merged$qseqid))) %>%
  filter(n <= 2)

# read in genome
genome_seq <- Biostrings::readDNAStringSet(filepath = paste0("seq/", opt$genome_name))
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# read in repeat library
consensus_seq <- Biostrings::readDNAStringSet(filepath = paste0("data/", opt$repeatmodeler))
names(consensus_seq) <- sub(" .*", "", names(consensus_seq))

# create preliminary multiple sequence fastas for self blast
for(i in 1:length(blast_out_tbl$qseqid)){

  align_ranges <- blast_out_merged[blast_out_merged$qseqid == blast_out_tbl$qseqid[i]]
  align_ranges <- GenomicRanges::reduce(align_ranges, ignore.strand=T)
  align_seq <- BSgenome::getSeq(genome_seq, align_ranges)
  names(align_seq) <- base::paste0(GenomeInfoDb::seqnames(align_ranges), ":", IRanges::ranges(align_ranges))
  
  align_seq <- c(consensus_seq[names(consensus_seq) == blast_out_tbl$qseqid[i]], align_seq)
  
  Biostrings::writeXStringSet(x = align_seq, filepath = paste0("out/initial_seq/", opt$genome_name, "_",
                                                               sub("#.*", "", blast_out_tbl$qseqid[i]), ".fasta"))

}