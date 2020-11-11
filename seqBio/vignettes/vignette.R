## ---- eval = FALSE------------------------------------------------------------
#  #Comprobación de los inputs:
#  stopifnot(class(seq1) == "character", class(seq2) == "character",
#            seq_type %in% c("protein", "dna"), seq_align  %in% c("local", "global"),
#            class(gap) == "numeric", class(N) == "numeric", class(shuff) == "numeric",
#            length(gap) == 2, length(N) == 1, shuff %in% c(1, 2))

## ---- eval = FALSE------------------------------------------------------------
#  s1 = readDNAStringSet(seq1, "fasta")
#  s2 = readDNAStringSet(seq2, "fasta")

## ---- eval = FALSE------------------------------------------------------------
#  S <- pairwiseAlignment(s1, s2, type = seq_align, substitutionMatrix = mat,
#                         gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)

## ----include=FALSE------------------------------------------------------------
#hello

