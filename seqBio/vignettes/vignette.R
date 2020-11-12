## ---- eval = FALSE------------------------------------------------------------
#  #Comprobación de los inputs:
#  stopifnot(class(seq1) == "character", class(seq2) == "character",
#            seq_type %in% c("protein", "dna"), seq_align  %in% c("local", "global"),
#            class(gap) == "numeric", class(N) == "numeric", class(shuff) == "numeric",
#            length(gap) == 2, length(N) == 1, shuff %in% c(1, 2))

## ---- eval = FALSE------------------------------------------------------------
#  if(seq_type == "dna"){
#    s1 = readDNAStringSet(seq1, "fasta")
#    s2 = readDNAStringSet(seq2, "fasta")
#  } else{
#    s1 = readAAStringSet(seq1, "fasta")
#    s2 = readAAStringSet(seq2, "fasta")
#  }

## ---- eval = FALSE------------------------------------------------------------
#  S <- pairwiseAlignment(s1, s2, type = seq_align, substitutionMatrix = mat,
#                         gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)

## ---- include=FALSE-----------------------------------------------------------
library(seqBio)

## -----------------------------------------------------------------------------
generateSeqsWithMultinomialModel("ATTGGCACT",5)

## ---- eval=FALSE--------------------------------------------------------------
#  if(shuff == 1){
#    randomseqs <- generateSeqsWithMultinomialModel(as.character(s1[[1]]), N)
#    randomscores <- double(N)
#  
#    for (i in 1:N){
#      score <- pairwiseAlignment(randomseqs[i], s2, type = seq_align, substitutionMatrix = mat,
#                                 gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)
#      randomscores[i] <- score
#    }
#  
#    } else {
#      randomseqs <- generateSeqsWithMultinomialModel(as.character(s2[[1]]), N)
#      randomscores <- double(N)
#  
#      for (i in 1:N){
#        score <- pairwiseAlignment(s1, randomseqs[i], type = seq_align, substitutionMatrix = mat,
#                                   gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)
#        randomscores[i] <- score
#      }
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  x <- mean(randomscores)
#  s <- sd(randomscores)

## ---- eval=FALSE--------------------------------------------------------------
#  lambda <- 1.2825 / s
#  u <- x - (0.45*s)

## ---- eval=FALSE--------------------------------------------------------------
#  m <- width(s1)
#  n <- width(s2)
#  K <- exp( (lambda*u) ) / (m*n)

## ---- eval=FALSE--------------------------------------------------------------
#  S_prima <- (lambda * S) - log( (K*m*n) )

## ---- eval = FALSE------------------------------------------------------------
#  prob <- sum(randomscores >= S) / N

## -----------------------------------------------------------------------------
library(seqBio)
library(Biostrings)
seq1 <- "~/GitHub/Bioinformatics_Task3/seqBio/data/gi32141095_N_0.fa"
seq2 <- "~/GitHub/Bioinformatics_Task3/seqBio/data/gi32141095_N_1.fa"
seq_type ="dna"
seq_align = "global"
mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
gap = c(-4, -1)
N = 200
shuff = 2

## -----------------------------------------------------------------------------
seq_alignment(seq1, seq2, seq_type, seq_align, mat, gap, N, shuff)

## -----------------------------------------------------------------------------
seq1 <- "~/GitHub/Bioinformatics_Task3/seqBio/data/P0DP27.fa"
seq2 <- "~/GitHub/Bioinformatics_Task3/seqBio/data/Q9N1R0.fa"
seq_type ="protein"
seq_align = "local"
mat = data("BLOSUM62")
gap = c(-12, -2)
N = 1000
shuff = 1

## -----------------------------------------------------------------------------
seq_alignment(seq1, seq2, seq_type, seq_align, mat, gap, N, shuff)

