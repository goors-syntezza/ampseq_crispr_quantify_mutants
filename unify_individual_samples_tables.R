#!/usr/bin/Rscript
fns <- Sys.glob("*/grna_mutation_rates_and_counts_per_amplicon.tsv")
fns <- fns[fns != 'wt/grna_mutation_rates_and_counts_per_amplicon.tsv']
my_df <- read.table(fns[1], sep = "\t", head = T, quote = "", stringsAsFactors = F)
for (i in 2 : length(fns)) {curr_df <- read.table(fns[i], sep = "\t", head = T, quote = "", stringsAsFactors = F);my_df <- rbind(my_df, curr_df)}
write.table(my_df, file = 'grna_mutation_rates_and_counts_per_amplicon_all_samples_all.tsv', sep = "\t", quote = F, col.names = T, row.names = F)

fns <- Sys.glob("*/ampliconDistribution_w_gRNA_region_and_mutation_status.tsv")
fns <- fns[fns != 'wt/ampliconDistribution_w_gRNA_region_and_mutation_status.tsv']
my_df <- read.table(fns[1], sep = "\t", head = T, quote = "", stringsAsFactors = F)
for (i in 2 : length(fns)) {curr_df <- read.table(fns[i], sep = "\t", head = T, quote = "", stringsAsFactors = F);my_df <- rbind(my_df, curr_df)}
my_df$seq_id <- paste0('seq',as.numeric(factor(my_df[, 'amplicon_seq'])))
write.table(my_df, file = 'ampliconDistribution_w_gRNA_region_and_mutation_status_all.tsv', sep = "\t", quote = F, col.names = T, row.names = F)
