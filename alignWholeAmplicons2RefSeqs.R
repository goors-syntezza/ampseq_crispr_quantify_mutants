#!/usr/bin/Rscript
library(Biostrings)
min_count = 10
alignSeq2RefAndGetAlignedRegionWithoutInsertions <- function(query, subject) {

 alignResult <- pairwiseAlignment(query, subject, type = "global")

 subject_aligned <- as.character(subject(alignResult))
 query_aligned <- as.character(pattern(alignResult))
 #print(subject_aligned)
 #print(query_aligned)
 alignment_length1 <- nchar(subject_aligned)
 alignment_length2 <- nchar(query_aligned)
 if (alignment_length1 != alignment_length2) stop("Error: expected aligned query and subject to be of same length!
")
 aligned_query_no_inserts <- ''
 for (i in 1 : alignment_length1) {
  currCharInAlignedSubject <- substring(subject_aligned, i, i)
  if (currCharInAlignedSubject == '-') {next} #print('#');
  currCharInAlignedQuery <- substring(query_aligned, i, i)
  aligned_query_no_inserts <- paste0(aligned_query_no_inserts, currCharInAlignedQuery)
  #print('*')
  #print(currCharInAlignedQuery)
  }
 print('**')
 if (nchar(aligned_query_no_inserts) > nchar(query)) {
  print('ERR');
  print(query)
  print(aligned_query_no_inserts)
  stop(1)
  
  }
 if (nchar(aligned_query_no_inserts) != nchar(query)) {
  print(alignResult@pattern@range@start)
  print(alignResult@pattern@range@width)
  startPos <- alignResult@pattern@range@start
  endPos <- startPos + alignResult@pattern@range@width
  print(endPos)
  print(nchar(query) - endPos)
  aligned_query_no_inserts <- paste0(paste0(rep('-', startPos -1)), aligned_query_no_inserts, paste0(rep('-', max(0, nchar(query) - endPos))))
  }
 
 return(aligned_query_no_inserts)
 
 }

my_amplicons_table <- read.table('amplicon_abundance_table_w_primer_ids.tsv', sep = "\t", quote = "", stringsAsFactors = F, head = F)
colnames(my_amplicons_table) <- c('Count', 'Seq', 'amplicon_id')
my_amplicons_table <- my_amplicons_table[my_amplicons_table$Count >= min_count, ]
print(dim(my_amplicons_table))
ref_amplicons <- read.table('../00.Metadata/amplicon_sequences.tsv', sep = "\t", quote = "", stringsAsFactors = F, head = T)

for (amplicon_id in ref_amplicons$amplicon_id) {
 my_amplicons_table_curr <- my_amplicons_table[my_amplicons_table$amplicon_id == amplicon_id, ]
 my_ref_seq <- ref_amplicons[ref_amplicons$amplicon_id == amplicon_id, 'seq']
 my_amplicons_table_curr$alignment <- sapply(my_amplicons_table_curr$Seq, function(x) alignSeq2RefAndGetAlignedRegionWithoutInsertions(query = x, subject = my_ref_seq))
 write.table(my_amplicons_table_curr, file = paste0(amplicon_id, '_whole_amplicon_alignment.tsv'), sep = "\t", quote = F, col.names = T, row.names = F)
 }