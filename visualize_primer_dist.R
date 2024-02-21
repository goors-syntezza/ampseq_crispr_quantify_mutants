#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {print("Error: please provide sample_id parameter!");stop(1)}
sample_id = args[1]
amplicon_abundance_table <- read.table(file = paste0('amplicon_abundance_table_w_primer_ids.tsv'), head = F, stringsAsFactors = F, sep = "\t", quote = "")
amplicon_abundance_table_agg <- aggregate(amplicon_abundance_table[, 1], by = list(amplicon_abundance_table[, 3]), FUN = sum)
colnames(amplicon_abundance_table_agg) <- c('primer_id', 'seq_count')
amplicon_abundance_table_agg$primer_id <- as.character(amplicon_abundance_table_agg$primer_id)
write.table(amplicon_abundance_table_agg, file = 'primers_dist_table.tsv', sep = "\t", quote = F, col.names = T, row.names = F)
png(paste0('primers_dist.png'))
pie(labels = paste0(amplicon_abundance_table_agg$primer_id, '\n', amplicon_abundance_table_agg$seq_count), x = amplicon_abundance_table_agg$seq_count, main = paste0(sample_id, ': amplicons read distribution'))
dev.off()
