#!/usr/bin/Rscript

my_tab_fn <- 'ampliconDistribution_w_gRNA_region_and_mutation_status.tsv'
my_tab <- read.table(my_tab_fn, sep = "\t", head = T, quote = "", stringsAsFactors = F)
my_tab <- my_tab[, c('amplicon_id', 'gRNA_id', 'gRNA_region_seq')]
my_tab <- my_tab[!duplicated(my_tab), ]

write.table(my_tab, file = 'wild_type_ref_grna_region_seqs.tsv', sep = "\t", quote = F, col.names = T, row.names = F)