#!/usr/bin/Rscript
library(Biostrings)
alignSeq2RefAndGetAlignedRegionWithoutInsertions <- function(query, subject) {
 
 alignResult <- pairwiseAlignment(query, subject, type = "global")
 
 subject_aligned <- as.character(subject(alignResult))
 query_aligned <- as.character(pattern(alignResult))
 print(subject_aligned)
 print(query_aligned)
 alignment_length1 <- nchar(subject_aligned)
 alignment_length2 <- nchar(query_aligned)
 if (alignment_length1 != alignment_length2) stop("Error: expected aligned query and subject to be of same length!")
 aligned_query_no_inserts <- ''
 for (i in 1 : alignment_length1) {
  currCharInAlignedSubject <- substring(subject_aligned, i, i)
  if (currCharInAlignedSubject == '-') {print('#');next}
  currCharInAlignedQuery <- substring(query_aligned, i, i)
  aligned_query_no_inserts <- paste0(aligned_query_no_inserts, currCharInAlignedQuery) 
  print('*')
  print(currCharInAlignedQuery)
  }
 return(aligned_query_no_inserts)
 }
 
print(alignSeq2RefAndGetAlignedRegionWithoutInsertions("AACTGACCC", "ACTGAACCC"))
print(alignSeq2RefAndGetAlignedRegionWithoutInsertions("AACTGAAAAAAAAAACTC", "ACTGAACCC"))
