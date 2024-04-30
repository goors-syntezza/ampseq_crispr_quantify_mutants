#!/usr/bin/Rscript
sink(stdout(), type="message")
library(ggplot2)
library(forcats)
library(rlist)

args = commandArgs(trailingOnly=TRUE)
metadata_dir <- '../00.Metadata/'
if (length(args) != 2) {print("Error: please provide sample_id, hasWTSample {True, False}!");stop(1)}
sample_id = args[1]
hasWTSample = args[2]
if ((sample_id != 'wt') && (hasWTSample == 'True')) wt_grna_region_seqs_df <- read.table('../wt/wild_type_ref_grna_region_seqs.tsv', sep = "\t", quote = "", head = T, stringsAsFactors = F)

is_grna_in_ref_seqs_provided <- function(grna_region_seq, amplicon_id) {
 exists_in_ref <- F
 my_wt_seqs <- unname(unlist(wt_grnas_per_exon_tab[wt_grnas_per_exon_tab$amplicon_id == amplicon_id, 'seq']))
 for (my_wt_ref_seq in my_wt_seqs) if (grepl(grna_region_seq, my_wt_ref_seq)) {exists_in_ref <- T;break}
 return(exists_in_ref)
 }

is_grna_in_wt <- function(grna_region_seq, amplicon_id, grna_id) {
 #print(grna_region_seq)
 my_wt_seqs <- unname(unlist(wt_grna_region_seqs_df[wt_grna_region_seqs_df$amplicon_id == amplicon_id & wt_grna_region_seqs_df$gRNA_id == grna_id, 'gRNA_region_seq']))
 #print(my_wt_seqs)
 is_in_wt_ref <- unname(grna_region_seq) %in% my_wt_seqs
 #print(is_in_wt_ref)
 return(is_in_wt_ref)
 }



seqs_dist_fn <- paste0('ampliconDistribution_w_gRNA_region.tsv')
seqs_dist_tab <- read.table(seqs_dist_fn, sep = "\t", quote = "", stringsAsFactors = F, head = T)

wt_grnas_per_exon_fn <- paste0(metadata_dir, 'amplicon_sequences.tsv')
wt_grnas_per_exon_tab <- read.table(wt_grnas_per_exon_fn, sep = "\t", quote = "", stringsAsFactors = F, head = T)

is_mutated <- apply(seqs_dist_tab, 1, function(x) {
 gRNA_id <- x[1]
 gRNA_region_seq <- x[8]
 amplicon_id <- x[4]
 if ((sample_id == 'wt') || (hasWTSample == 'False')) ret_val <- !is_grna_in_ref_seqs_provided(grna_region_seq = gRNA_region_seq, amplicon_id = amplicon_id)
 else ret_val <- !is_grna_in_wt(grna_region_seq = gRNA_region_seq, amplicon_id = amplicon_id, grna_id = gRNA_id)
 return(ret_val)
 })

seqs_dist_tab$is_mutated <- is_mutated
write.table(seqs_dist_tab, file = 'ampliconDistribution_w_gRNA_region_and_mutation_status.tsv', sep = "\t", quote = F, col.names = T, row.names = F)
#colnames(my_results_as_df) <- c('sample_id', 'exon_id', 'wt_reads_count', 'mut_reads_count', 'percent_mut')
#print(dim(my_results_as_df))
#write.table(my_results_as_df, file = 'mutation_rate_per_sample_and_exon_based_on_intact_gRNA_plus_flank5_presence.tsv', sep = "\t", col.names = T, row.names = F, quote = F)
warnings()
