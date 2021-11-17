#!/usr/bin/env Rscript

# suppressPackageStartupMessages({
#   library(optparse)
# })
# 
# 
# # parse input variables
# option_list = list(
#   make_option(c("-g", "--genome_name"), type="character", default=NULL, 
#               help="genome name", metavar="character"),
#   make_option(c("-s", "--species_name"), type="character", default=NULL,
#               help="genome name", metavar="character")
# )
# 
# opt_parser = OptionParser(option_list=option_list)
# opt = parse_args(opt_parser)
# 
# if (is.null(opt$genome_name)) {
#   stop("Genome name is needed")
# } else {
#   # set genome names
#   genome_name <- opt$genome_name
# }
# 
# if (is.null(opt$species_name)) {
#   species_name <- opt$genome_name
# } else {
#   species_name <- opt$species_name
# }

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

flank5 = 2500
flank3 = 2500
genome_name = "myzusPersicae_bassLabG006.fasta"

queries <- read_tsv(paste0("out/", genome_name, "_self_queries.txt"), col_names = "qseqid")
i=857
for (i in 1:nrow(queries)) {
  
  in_seq <- Biostrings::readDNAStringSet(filepath = paste0("out/initial_seq/", queries$qseqid[i]))
  
  self_out <- read_tsv(file = paste0("out/self_search/", queries$qseqid[i], ".out"),
                       col_names = c("seqnames", "sseqid", "pident", "length", "qstart", "qend",
                                     "qlen", "sstart", "send", "slen", "evalue", "bitscore")) %>%
    filter(sseqid != seqnames) %>%
    filter(!grepl("rnd-._family-.*#", seqnames)) %>%
    filter(sub("#.*", "", seqnames) != queries$qseqid[i]) %>%
    filter(length > 0.5 * BiocGenerics::width(in_seq[1]))
  
  self_ranges <- self_out %>%
    group_by(seqnames) %>%
    mutate(start = min(qstart),
           end = max(qend)) %>%
    dplyr::ungroup() %>%
    dplyr::select(seqnames, start, end) %>%
    base::unique() %>%
    dplyr::arrange(seqnames, start, end) %>%
    plyranges::as_granges()
    
  out_seq <- getSeq(in_seq, self_ranges)
  names(out_seq) <- seqnames(self_ranges)
  
  out_seq <- c(in_seq[1], out_seq)
  
  writeXStringSet(out_seq, paste0("out/to_align/", queries$qseqid[i]))

}
