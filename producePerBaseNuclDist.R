#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {print("Error: please provide sample_id parameter!");stop(1)}
sample_id = args[1]


gRNAalignmentsAndCount <- read.table('ampliconDistribution_w_gRNA_region_and_mutation_status.tsv', head = T, sep = "\t", quote = "", stringsAsFactors = F)

AmpGRNaCombs <- gRNAalignmentsAndCount[, c('amplicon_id', 'gRNA_id')]
AmpGRNaCombs <- AmpGRNaCombs[!duplicated(AmpGRNaCombs), ]

for (i in 1 : nrow(AmpGRNaCombs)) {
 curr_amp_id <- AmpGRNaCombs[i, 'amplicon_id']
 curr_grna_id <- AmpGRNaCombs[i, 'gRNA_id']
 
 curr_data <- gRNAalignmentsAndCount[gRNAalignmentsAndCount$amplicon_id == curr_amp_id & gRNAalignmentsAndCount$gRNA_id == curr_grna_id, c('count', 'gRNA_region_seq_globalAlign')]
 
 grna_len <- nchar(curr_data[1 , 'gRNA_region_seq_globalAlign'])
 
 basesCountMat <- data.frame(matrix(0, nrow = grna_len, ncol = 7))
 colnames(basesCountMat) <- c('Major', 'A', 'C', 'T', 'G', 'N', 'Del')
 
 for (i in 1 : nrow(curr_data)) for (j in 1: grna_len)
  {
   my_count <- curr_data[i, 'count']
   my_base <- substr(curr_data[i, 'gRNA_region_seq_globalAlign'], j, j)
   #print(my_base)
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
 write.table(basesCountMat, file = paste0('gRNAbaseCountPerPosition_', curr_amp_id ,'_', curr_grna_id, '_', sample_id, '.tsv'), sep = "\t", quote = F, col.names = T, row.names = F)
 write.table(basesCountMatPercent, file = paste0('gRNAbasePercentPerPosition_', curr_amp_id ,'_', curr_grna_id, '_', sample_id, '.tsv'), sep = "\t", quote = F, col.names = T, row.names = F)
 }