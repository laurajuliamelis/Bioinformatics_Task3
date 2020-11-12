#' Sequence alignment and statistical significance.
#'
#' La funcion \code{seq_alignment} muestra la significacion estadistica de un alineamiento
#' por parejas.
#' 
#' @param seq1 Secuencia 1 en formato Fasta.
#' @param seq2 Secuencia 2 en formato Fasta.
#' @param seq_type Tipo de secuencia: "protein" o "dna"
#' @param seq_align Tipo de alineamiento: "local" o "global".
#' @param mat Matriz de substitución: PAMn, BLOSUMn, … Ver \code{data(package="Biostrings")} para ver todas las opciones disponibles.
#' @param gap Puntuación Gap: Open penalty, extended penalty.
#' @param N Numero de replicas.
#' @param shuff Secuencia donde se hace shuffling: 1 o 2.
#'
#' @return  \code{seq_alignment} devuelve un resultado gráfico (histograma de las puntiaciones estimadas) y un resultado numérico.
#'
#' @examples
#' 
#' seq_alignment(seq1 = "P0DP27.fa", seq2 = "Q9N1R0.fa", seq_type ="protein",
#'               seq_align = "local", mat = data("BLOSUM62"), gap = c(-12, -2), N = 100, shuff = 1)
#'               
#' mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
#' seq_alignment(seq1 = "gi32141095_N_0.fa", seq2 = "gi32141095_N_1.fa", seq_type ="dna",
#'               seq_align = "global", mat = mat, gap = c(-8, -1), N = 2000, shuff)
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
  
  if(seq_type == "dna"){
    s1 = readDNAStringSet(seq1, "fasta")
    s2 = readDNAStringSet(seq2, "fasta")
    name_mat = "DNAScoringMat"
  } else{
    s1 = readAAStringSet(seq1, "fasta")
    s2 = readAAStringSet(seq2, "fasta")
    name_mat = mat
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
  ###      puntuación superior a la putuación original
  prob <- sum(randomscores >= S) / N
  
  # 6. RESULTADO GRÁFICO
  require("extRemes", quietly = T)
  require("gplots", quietly = T) 

  fit0 <- fevd(randomscores, type = "Gumbel")
  histo <- hist(randomscores, breaks = 20, plot= FALSE)
  y <- histo$density[max(which(histo$mids <= S))]
  
  texto <- c(sub(".*/([^.]+)\\..*", "\\1", seq1), sub(".*/([^.]+)\\..*", "\\1", seq2), seq_type, seq_align,
             name_mat, as.character(paste(gap[1], "y", gap[2])), N, shuff,
             round(lambda,4), round(u,4), round(K,4), S, round(S_prima,4), round(prob,4))
  texto <- as.data.frame(texto)
  rownames(texto) <- c("Nombre secuencia 1:", "Nombre secuencia 2:", "Tipo de alineamiento:", 
                       "Tipo de secuencia;", "Matriz de substitución",  "Puntuación Gap:",
                       "Número de réplicas:", "Secuencia shuffling:", "Parámetro lambda:", 
                       "Parámetro u:", "Constante K:", "Score original:", "Score estandarizado (S'):", 
                       "P(x >= S)")
  colnames(texto) <- " "
  
#  if(.Platform$OS.type == "unix") {
#    quartz(width = 10, height = 7)
#  } else {
#    windows(width = 10, height = 7)
#  }
  
  par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
  plot(fit0, "hist", ylim=c(min(density(randomscores)$y),max(density(randomscores)$y)+0.05), 
       col="lightblue", lty= 1, lwd = 1, main = " ")
  text(x = S, y = 0.003, labels = paste0("S=", S), cex = 0.75, font=2, col ="red")
  text(x = S, y = 0.001, labels = sprintf('\u2193'), cex = 1, font=2, col ="red")
  
  textplot(texto, valign = "top", halign = "right", cex =0.80)
  title(paste("Histograma de los N =", N,"scores"), outer = TRUE, cex.main = 1.5)
  
  # 7. REUSLTADO NUMÉRICO
  txt1 <- "Resumen numérico de las puntuaciones obtenidas al hacer shuffling:"
  cat(txt1, fill = TRUE)
  cat(rep("=", nchar(txt1)), sep = "", fill = TRUE)
  print(summary(randomscores))
  
  cat("\n\n")
  txt2 <- "Otros resultados numéricos:"
  cat(txt2, fill = TRUE)
  cat(rep("=", nchar(txt2)), sep = "", fill = TRUE)
  nums <- c(round(lambda,2) , round(u,2), round(K,2), round(S,2), round(S_prima,2), round(prob,2))
  nums <- as.data.frame(nums)
  rownames(nums) <- c("lambda", "u", "K", "S", "S'", "P(x >= S)")
  colnames(nums) <- NULL
  print(nums)
  cat("\n\n")
}