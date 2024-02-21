#!/usr/bin/Rscript
library(reshape2)

addMissingZeroCounts <- function(my_tab) {
 orig_nrow <- nrow(my_tab)
 for (i in 1 : orig_nrow) {
  sample_id <- my_tab[i, 'sample_id']
  amplicon_id <- my_tab[i, 'amplicon_id']
  gRNA_id <- my_tab[i, 'gRNA_id']
  comp_mut_status <- as.logical(my_tab[i, 'is_mutated'])
  comp_mut_status <- !comp_mut_status
  if (nrow(my_tab[my_tab$sample_id == sample_id & my_tab$amplicon_id == amplicon_id & my_tab$gRNA_id == gRNA_id & my_tab$gRNA_id == gRNA_id & my_tab$is_mutated == comp_mut_status, ]) == 0) my_tab <- rbind(my_tab, c(sample_id, amplicon_id, gRNA_id, comp_mut_status, 0))
  }
 return(my_tab)
 }

my_tab <- read.table('ampliconDistribution_w_gRNA_region_and_mutation_status.tsv', sep = "\t", head = T, quote = "", stringsAsFactors = F)
my_tab_agg <- aggregate(count ~ sample_id + amplicon_id + gRNA_id + is_mutated, data = my_tab, FUN = sum, na.rm = TRUE)
my_tab_agg <- addMissingZeroCounts(my_tab_agg)
my_tab_agg_wide <- dcast(my_tab_agg, sample_id + amplicon_id + gRNA_id ~ is_mutated, value.var="count")
colnames(my_tab_agg_wide)[4:5] <- c('wt_reads_count', 'mut_reads_count')
my_tab_agg_wide$mut_reads_count <- as.numeric(my_tab_agg_wide$mut_reads_count)
my_tab_agg_wide$wt_reads_count  <- as.numeric(my_tab_agg_wide$wt_reads_count)
my_tab_agg_wide$percent_mut <- my_tab_agg_wide$mut_reads_count / (my_tab_agg_wide$mut_reads_count + my_tab_agg_wide$wt_reads_count) * 100
write.table(my_tab_agg_wide, file = 'grna_mutation_rates_and_counts_per_amplicon.tsv', sep = "\t", col.names = T, quote = F, row.names = F)