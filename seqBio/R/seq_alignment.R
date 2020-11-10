#' Sequence alignment and statistical significance.
#'
#' La función \code{seq_alignment} muestra la significación estadística de un alineamiento
#' por parejas.
#' 
#' @param seq1 Secuencia 1 en formato Fasta.
#' @param seq2 Secuencia 2 en formato Fasta.
#' @param seq_type Tipo de secuencia: Proteina o ADN.
#' @param seq_align Tipo de alineamiento: "local" o "global".
#' @param mat Matriz de substitución: PAMn, BLOSUMn, … Ver \code{data(package="Biostrings")} para ver todas las opciones disponibles.
#' @param gap Puntuación Gap: Open penalty, extended penalty.
#' @param N Número de replicas.
#' @param shuff Secuencia donde se hace shuffling: 1 o 2.
#'
#' @return  \code{seq_alignment} devuelve un resultado gráfico (histograma de las puntiaciones estimadas) y un resultado numérico.
#'
#' @examples
#' 
#' seq_alignment(seq1 = "P0DP27.fa", seq2 = "Q9N1R0.fa", seq_type ="protein",
#'               seq_align = "local", mat = data("BLOSUM62"), gap = c(-12, -2), N = 100, shuff = 1)
#' seq_alignment(seq1 = "gi32121095_N_0.fa", seq2 = "gi32121095_N_1.fa", seq_type ="dna",
#'               seq_align = "global", mat = data("BLOSUM6100"), gap = c(-8, -1), N = 2000, shuff)
#' seq_alignment(seq1 = "P69905.fa", seq2 = "P01942.fa", seq_type ="protein",
#'               seq_align = "global", mat = data("BLOSUM62"), gap = c(10, 0.5), N =1000, shuff=1)
#' 
#' @references \url{https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter4.html}
#'
#' @export
#' 

seq_alignment <- function (seq1, seq2, seq_type = c("protein", "dna"), seq_align = c("local", "global"), mat, gap, N, shuff = c(1, 2)){
  # 1. Carga del paquete Biostrings:
  if(!require('Biostrings',quietly = T)){BiocManager::install("Biostrings")}
  require('Biostrings', quietly = T)
  
  # 2. Comprobación de los inputs:
  stopifnot(class(seq1) == "character", class(seq2) == "character",
            seq_type %in% c("protein", "dna"), seq_align  %in% c("local", "global"),
            class(gap) == "numeric", class(N) == "numeric", class(shuff) == "numeric",
            length(gap) == 2, length(N) == 1, shuff %in% c(1, 2))
  
  # 3. Lectura de las secuencias:
  s1 = readDNAStringSet(seq1, "fasta")
  s2 = readDNAStringSet(seq2, "fasta")
  
  # 4. Creación de una funcion de haga N shufflings:
  
  generateSeqsWithMultinomialModel <- function(inputsequence, X){
    
    # Change the input sequence into a vector of letters
    require("seqinr", quietly = T) # This function requires the SeqinR package.
    inputsequencevector <- s2c(inputsequence)
    
    # Find the frequencies of the letters in the input sequence "inputsequencevector":
    mylength <- length(inputsequencevector)
    mytable <- table(inputsequencevector)
    
    # Find the names of the letters in the sequence
    letters <- rownames(mytable)
    numletters <- length(letters)
    probabilities <- numeric() # Make a vector to store the probabilities of letters
    for (i in 1:numletters){
      letter <- letters[i]
      count <- mytable[[i]]
      probabilities[i] <- count/mylength
    }
    
    # Make X random sequences using the multinomial model with probabilities "probabilities"
    seqs <- numeric(X)
    for (j in 1:X){
      seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
      seq <- c2s(seq)
      seqs[j] <- seq
    }
    
    # Return the vector of random sequences
    return(seqs)
  }
  
  # 4. Alineamiento original (óptimo):
  S <- pairwiseAlignment(s1, s2, type = seq_align, substitutionMatrix = mat, 
                             gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)  
  
  # 5. Shuffling:
  
  if(shuff == 1){
    randomseqs <- generateSeqsWithMultinomialModel(as.character(s1[[1]]), N)
    randomscores <- double(N)
    for (i in 1:N){
      score <- pairwiseAlignment(randomseqs[i], s2, type = seq_align, substitutionMatrix = mat, 
                                 gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)
      randomscores[i] <- score
    }
    
  } else{
    randomseqs <- generateSeqsWithMultinomialModel(as.character(s2[[1]]), N)
    randomscores <- double(N)
    
    for (i in 1:N){
      score <- pairwiseAlignment(s1, randomseqs[i], type = seq_align, substitutionMatrix = mat, 
                                 gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)
      randomscores[i] <- score
    }
    
  }
  # 6. Cálculos:
  
  ### 6.1. Media y desv estándar de las puntiaciones.
  x <- mean(randomscores)
  s <- sd(randomscores)
  
  ### 6.2. Estimación de los parámetros de la distribución Gumbel.
  lambda <- 1.2825 / s
  u <- x - (0.45*s)
  
  ### 6.3. Estimación de K
  m <- width(s1)
  n <- width(s2)
  K <- exp( (lambda*u) ) / (m*n)
  
  ### 6.4. Estandarización de la puntuación S
  S_prima <- (lambda * S) - log( (K*m*n) )
  
  ### 6.5. Probabilidad que un alineamiento de secuencias al azar tenga una
  ###      puntuación superior a la obtenida por azar
  prob <- 1 - exp( exp(-S_prima)) 
  
  # 6. RESULTADO GRÁFICO
  require("reliaR", quietly = T)
  require("gplots", quietly = T) 
  
  dens <- dgumbel(randomscores, mu = u, sigma = lambda, log = FALSE)
  gens <- lambda * exp((-lambda)*(randomscores-u) - exp((-lambda)*(randomscores-u)))
  histo <- hist(randomscores, col="lightblue", breaks = 20, freq = FALSE, plot= FALSE)
  y <- histo$density[max(which(histo$mids <= S))]
  
  texto <- c(gsub('.{3}$', '', seq1), gsub('.{3}$', '', seq2), seq_type, seq_align, mat, as.character(paste(gap[1], "y", gap[2])), N, shuff,
             round(lambda,4), round(u,4), round(K,4), S, round(S_prima,4), round(prob,4))
  texto <- as.data.frame(texto)
  rownames(texto) <- c("Nombre secuencia 1:", "Nombre secuencia 2:", "Tipo de alineamiento:", 
                       "Tipo de secuencia;", "Matriz de substitución",  "Puntuación Gap:",
                       "Número de réplicas;", "Secuencia shuffling:", "Parámetro lambda:", 
                       "Parámetro u:", "Constante K:", "Score original:", "Score estandarizado (S'):", 
                       "P(S')")
  colnames(texto) <- " "
  
  if(.Platform$OS.type == "unix") {
    quartz(width = 8, height = 7)
  } else {
    windows(width = 10)
  }
  
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
  hist(randomscores, col="lightblue", freq = FALSE, xaxt='n', main = NULL)
  axis(1, at = seq(range(randomscores)[1], range(randomscores)[2], by=10), las=1)
  lines(dens, lty= 1, lwd = 1, col= "red")
  text(x = S, y = (y+0.0115), labels = paste0("S=", S), cex = 0.75, font=2)
  text(x = S, y = (y + 0.005), labels = sprintf('\u2193'), cex = 1, font=2)
  
  textplot(texto, valign = "top", halign = "right", cex =0.80)
  title(paste("Histograma de los N =", N,"scores"), outer = TRUE, cex.main = 1.5)
  
  # 7. REUSLTADO NUMÉRICO
  cat("\n\n")
  txt1 <- "Resumen numérico de las puntuaciones obtenidas al hacer shuffling:"
  cat(txt1, fill = TRUE)
  cat(rep("=", nchar(txt1)), sep = "", fill = TRUE)
  print(summary(randomscores))
  
  cat("\n\n")
  txt2 <- "Otros resultados numéricos:"
  cat(txt2, fill = TRUE)
  cat(rep("=", nchar(txt2)), sep = "", fill = TRUE)
  nums <- c(lambda, u, K, S, S_prima, prob)
  nums <- as.data.frame(nums)
  rownames(nums) <- c("lambda", "u", "K", "S", "S'", "P(S' >= x)")
  colnames(nums) <- NULL
  print(nums)
  cat("\n\n")
}