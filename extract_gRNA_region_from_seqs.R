#!/usr/bin/Rscript
library(dplyr)
library(ggplot2)
library(forcats)
library(rlist)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {print("Error: please provide sample_id, min_seq_depth and flank_len parameters!");stop(1)}
sample_id = args[1]
min_seq_depth = as.numeric(args[2])
flank_len = as.numeric(args[3])

alignSeqToPatternAndReturnAlignedRegion <- function(my_pattern, my_seq, flankSeqLength) {
 alignResult <- pairwiseAlignment(my_pattern, my_seq, type = 'overlap')
 
 indices <- as.data.frame(alignResult@subject@range)
 start <- indices[['start']] - flankSeqLength
 end <- indices[['end']] + flankSeqLength
 alignRegion <- substr(my_seq, start, end)
 return(alignRegion)
 }

seqs_dist_fn <- paste0(sample_id, '_amplicon_abundance_table_w_primer_ids.tsv')
grna_seqs_fn <- 'gRNA_sequences.tsv'
seqs_dist_tab <- read.table(seqs_dist_fn, sep = "\t", quote = "", stringsAsFactors = F, head = T)
colnames(seqs_dist_tab) <- c('count', 'seq', 'amplicon_id')
seqs_dist_tab <- seqs_dist_tab[seqs_dist_tab$count >= min_seq_depth, ]
seqs_dist_tab <- seqs_dist_tab[seqs_dist_tab$amplicon_id != 'Unidentified', ]

nrow(seqs_dist_tab)

grna_seqs_tab <- read.table(grna_seqs_fn, sep = "\t", head = T, quote = "", stringsAsFactors = F)
print(grna_seqs_tab)

seq_indices <- paste0(sample_id, '_ampseq_', 1 : nrow(seqs_dist_tab), sep = "")
rownames(seqs_dist_tab) <- seq_indices
amplicon_ids <- seqs_dist_tab[, 'amplicon_id']

alignmanets2make_table <- lapply(seq_indices, function(x) {
 seq_index <- x
 seq_amplicon_id <- seqs_dist_tab[x, 'amplicon_id']
 seq_grnas <- grna_seqs_tab[grna_seqs_tab$amplicon_id == seq_amplicon_id, c('gRNA_id', 'seq'), drop = F]
 colnames(seq_grnas) <- c('gRNA_id', 'gRNA_seq')
 seq_grnas$seq_id <- x
 seq_grnas$amplicon_id <- seq_amplicon_id
 seq_grnas$amplicon_seq <- seqs_dist_tab[seq_index, 'seq']
 seq_grnas$sample_id <- sample_id
 seq_grnas$count <- seqs_dist_tab[seq_index, 'count']
 return(seq_grnas)
 })

#alignmanets2make_table <- data.frame(alignmanets2make_table)
alignmanets2make_table <- bind_rows(alignmanets2make_table)
print(alignmanets2make_table)
print(dim(alignmanets2make_table))

my_gRNA_Aligned_Regions <- apply(alignmanets2make_table, 1 ,function(x) {
	amplicon_id <- unname(x[4])
	seq_id <- unname(x[3])
	curr_gRNA_seq <- unname(x[2])
	amplicon_seq <-  unname(x[5])
	alignedRegionSeq <- alignSeqToPatternAndReturnAlignedRegion(curr_gRNA_seq, amplicon_seq, flank_len)
	return(alignedRegionSeq)
	})

my_gRNA_Aligned_Regions <- unname(unlist(my_gRNA_Aligned_Regions))
print(my_gRNA_Aligned_Regions)
alignmanets2make_table$gRNA_region_seq <- my_gRNA_Aligned_Regions
write.table(alignmanets2make_table, paste0('ampliconDistribution_w_gRNA_region_seq_flank_', flank_len,'.csv'), sep = "\t", quote = F, col.names = T, row.names = F)

