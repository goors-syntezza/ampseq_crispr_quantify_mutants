#!/usr/bin/Rscript
library(reshape2)
library(dplyr)
library(scatterpie)
library(ggplot2)
library(tidyverse)
library(gtools)
library(forcats)

mut_rates_and_reads_counts_tab_fn <- 'grna_mutation_rates_and_counts_per_amplicon_all_samples_all.tsv'
mut_rates_and_reads_counts_tab <- read.table(mut_rates_and_reads_counts_tab_fn, sep = "\t", quote = "", head = T, stringsAsFactors = F)
mut_rates_and_reads_counts_tab$exon_id <- paste0(mut_rates_and_reads_counts_tab$amplicon_id, '_', mut_rates_and_reads_counts_tab$gRNA_id)
exon_pos_tab <- data.frame(exon_id = mut_rates_and_reads_counts_tab$exon_id, exon_pos = 1 : length(mut_rates_and_reads_counts_tab$exon_id))
 
mut_reads_sum_paer_sample <- aggregate(mut_rates_and_reads_counts_tab$mut_reads_count, by=list(sample_id=mut_rates_and_reads_counts_tab$sample_id), FUN=sum)
colnames(mut_reads_sum_paer_sample)[2] <- 'mut_rates_and_reads_counts_tab'
mutated_sample_ids <- mut_reads_sum_paer_sample[mut_reads_sum_paer_sample$mut_rates_and_reads_counts_tab > 0, 'sample_id']
print(mut_reads_sum_paer_sample)
print(mixedsort(mutated_sample_ids))
mut_rates_and_reads_counts_tab <- mut_rates_and_reads_counts_tab[mut_rates_and_reads_counts_tab$sample_id %in% mutated_sample_ids, ]
mut_rates_and_reads_counts_tab$total_reads_count <- mut_rates_and_reads_counts_tab$wt_reads_count + mut_rates_and_reads_counts_tab$mut_reads_count
mut_rates_and_reads_counts_tab$percent_unmut <- 100 - mut_rates_and_reads_counts_tab$percent_mut
mut_rates_and_reads_counts_tab_long <- melt(mut_rates_and_reads_counts_tab, id.vars = c('sample_id', 'exon_id', 'total_reads_count'), measure.vars=c('percent_unmut', 'percent_mut')) #measure.vars=c('total_reads_count', 'percent_mut')
colnames(mut_rates_and_reads_counts_tab_long) <- c('sample_id', 'exon_id', 'total_reads_count', 'Group', 'value')
mut_rates_and_reads_counts_tab_long$r <- log10(mut_rates_and_reads_counts_tab_long$total_reads_count) / 10
### old: mut_rates_and_reads_counts_tab_long$y <- as.numeric(sub('.', '', mut_rates_and_reads_counts_tab_long$sample_id))
mut_rates_and_reads_counts_tab_long$sample_id <- fct_inorder(as.character(mut_rates_and_reads_counts_tab_long$sample_id))
mut_rates_and_reads_counts_tab_long$y <- as.numeric(mut_rates_and_reads_counts_tab_long$sample_id)
print(mut_rates_and_reads_counts_tab_long$y)
loc_trans_table <- data.frame(from = unique(print(mut_rates_and_reads_counts_tab_long$y)), to = rank(unique(mut_rates_and_reads_counts_tab_long$y)))
mut_rates_and_reads_counts_tab_long$y <- loc_trans_table[match(mut_rates_and_reads_counts_tab_long$y, loc_trans_table$from), 'to']
print(loc_trans_table)

mut_rates_and_reads_counts_tab_long$x <- exon_pos_tab[match(mut_rates_and_reads_counts_tab_long$exon_id, exon_pos_tab$exon_id), 'exon_pos']
mut_rates_and_reads_counts_tab_long$x <- mut_rates_and_reads_counts_tab_long$x * 3.75

print(mut_rates_and_reads_counts_tab_long) #[1 : 5,]
print(unique(mut_rates_and_reads_counts_tab_long$exon_id))
print(unique(mut_rates_and_reads_counts_tab_long$x))
print(unique(mut_rates_and_reads_counts_tab_long$sample_id))
print(unique(mut_rates_and_reads_counts_tab_long$y))
print(mut_rates_and_reads_counts_tab_long)
print('*')

p <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, fill = Group, r = r), data = mut_rates_and_reads_counts_tab_long, cols = "Group", long_format = TRUE) +
  scale_x_continuous(breaks = unique(mut_rates_and_reads_counts_tab_long$x), labels = unique(mut_rates_and_reads_counts_tab_long$exon_id)) +
  scale_y_continuous(breaks = unique(mut_rates_and_reads_counts_tab_long$y), labels = unique(mut_rates_and_reads_counts_tab_long$sample_id)) +
  xlab('Exon ID') + ylab('Sample ID') + labs(fill = "Sequence\nstate") +
  theme(axis.text = element_text(color = 'black', size = 20), axis.title = element_text(color = 'black', size = 24), legend.text = element_text(color = 'black', size = 20), legend.title = element_text(color = 'black', size = 24)) +
  geom_scatterpie_legend(seq(min(mut_rates_and_reads_counts_tab_long$r), max(mut_rates_and_reads_counts_tab_long$r), length = 4),  x = 4 * length(unique(mut_rates_and_reads_counts_tab_long$exon_id)) + 2, y = length(unique(mut_rates_and_reads_counts_tab_long$sample_id)) - 2, labeller = function(x) 10^(x * 10)) +
  guides(scatterpie = guide_legend(override.aes = list(size=30)))

png('mutation_rate_and_reads_count_per_sample_and_exon_pie_lattice_based_on_intact_gRNA_region_only_mutated_samples.png', width = 1480, height = 1080)
print(p)
dev.off()
   
