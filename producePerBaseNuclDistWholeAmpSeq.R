#!/usr/bin/Rscript
min_count = 10
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {print("Error: please provide sample_id parameter!");stop(1)}
sample_id = args[1]
ref_amplicons <- read.table('../00.Metadata/amplicon_sequences.tsv', sep = "\t", quote = "", stringsAsFactors = F, head = T)

print(ref_amplicons)
for (amplicon_id in ref_amplicons$amplicon_id) {
 ref_seq <- ref_amplicons[ref_amplicons$amplicon_id == amplicon_id, 'seq']
 amplicon_len <- nchar(ref_seq)
 print(ref_seq)
 amplicon_seq_dist_and_alignment_fn <- paste0(amplicon_id ,'_whole_amplicon_alignment.tsv')
 my_amplicons_table <- read.table(amplicon_seq_dist_and_alignment_fn, sep = "\t", quote = "", stringsAsFactors = F, head = T)
 
 curr_data <- my_amplicons_table
 
 
 
 basesCountMat <- data.frame(matrix(0, nrow = amplicon_len, ncol = 8))
 colnames(basesCountMat) <- c('Ref', 'Major', 'A', 'C', 'T', 'G', 'N', 'Del')
 
 for (i in 1 : nrow(curr_data)) for (j in 1: amplicon_len)
  {
   my_count <- curr_data[i, 'Count']
   my_base <- substr(curr_data[i, 'alignment'], j, j)
   
   
   if (my_base == '-') my_base <- 'Del'
   
   basesCountMat[j, my_base] <- basesCountMat[j, my_base] + my_count
  }
 
 for (k in 1 : nrow(basesCountMat)) {
  max_base_index <-  which.max(basesCountMat[k, c('A', 'C', 'T', 'G', 'N', 'Del')])
  max_base_id <- c('A', 'C', 'T', 'G', 'N', 'Del')[max_base_index]
  max_base_value <- basesCountMat[k, max_base_id]
  MaxCount <- sum(basesCountMat[k, c('A', 'C', 'T', 'G', 'N', 'Del')] == max_base_value)
  if (MaxCount == 1) basesCountMat[k, 'Major'] <- max_base_id else basesCountMat[k, 'Major'] <- '-' 
  }
 
 basesCountMatPercent <- basesCountMat
 
 for (l in 1 : nrow(basesCountMatPercent)) {
  total_count <- sum(basesCountMatPercent[l , c('A', 'C', 'T', 'G', 'N', 'Del')])
  for (curr_base in c('A', 'C', 'T', 'G', 'N', 'Del')) basesCountMatPercent[l, curr_base] <- basesCountMatPercent[l, curr_base] / total_count * 100
  }
 
 #print(basesCountMat)
 #print(basesCountMatPercent)
 basesCountMat$Pos <- 1:nrow(print(basesCountMat))
 basesCountMatPercent$Pos <- 1:nrow(print(basesCountMatPercent))
 basesCountMat$Ref <- sapply(1 : length(ref_seq), function (i) substr(ref_seq, i, i))
 basesCountMatPercent$Ref <- sapply(1 : length(ref_seq), function (i) substr(ref_seq, i, i))
 write.table(basesCountMat, file = paste0('AmpliconSeq_', amplicon_id, '_CountPerPosition_', sample_id, '.tsv'), sep = "\t", quote = F, col.names = T, row.names = F)
 write.table(basesCountMatPercent, file = paste0('AmpliconSeq_', amplicon_id, '_PercentPerPosition_', sample_id, '.tsv'), sep = "\t", quote = F, col.names = T, row.names = F)
 }