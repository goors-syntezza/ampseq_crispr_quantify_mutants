#!/usr/bin/Rscript
sink(stdout(), type="message")
library(ggplot2)
library(forcats)
library(rlist)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {print("Error: please provide sample_id, flank_len parameters!");stop(1)}

is_grna_in_wt <- function(grna_region_seq, amplicon_id) {
 exists_in_ref <- F
 my_wt_seqs <- unname(unlist(wt_grnas_per_exon_tab[wt_grnas_per_exon_tab$amplicon_id == amplicon_id, 'seq']))
 for (my_wt_ref_seq in my_wt_seqs) if (grepl(grna_region_seq, my_wt_ref_seq)) {exists_in_ref <- T;break}
 return(exists_in_ref)
 }

sample_id = as.numeric(args[1])
flank_len = as.numeric(args[2])

seqs_dist_fn <- paste0('ampliconDistribution_w_gRNA_region_seq_flank_', flank_len, '.csv')
seqs_dist_tab <- read.table(seqs_dist_fn, sep = "\t", quote = "", stringsAsFactors = F, head = T)

wt_grnas_per_exon_fn <- 'amplicon_ref_seqs.tsv'
wt_grnas_per_exon_tab <- read.table(wt_grnas_per_exon_fn, sep = "\t", quote = "", stringsAsFactors = F, head = T)

is_mutated <- apply(seqs_dist_tab, 1, function(x) {
 gRNA_region_seq <- x[8]
 amplicon_id <- x[4]
 return(!is_grna_in_wt(grna_region_seq = gRNA_region_seq, amplicon_id = amplicon_id))
 })

seqs_dist_tab$is_mutated <- is_mutated
write.table(seqs_dist_tab, file = paste0('ampliconDistribution_w_gRNA_region_seq_flank_', flank_len, '_and_mutation_status.csv'), sep = "\t", quote = F, col.names = T, row.names = F)
#colnames(my_results_as_df) <- c('sample_id', 'exon_id', 'wt_reads_count', 'mut_reads_count', 'percent_mut')
#print(dim(my_results_as_df))
#write.table(my_results_as_df, file = 'mutation_rate_per_sample_and_exon_based_on_intact_gRNA_plus_flank5_presence.tsv', sep = "\t", col.names = T, row.names = F, quote = F)
warnings()
