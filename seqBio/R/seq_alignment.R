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
#' @return  \code{knapsack_dynamic} returns a list with two elements: the elements added to the knapsack and the maximum knapsack value.
#'
#' @examples
#' seq_alignment()
#' seq_alignment()
#' seq_alignment()
#' seq_alignment()
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
    require("seqinr") # This function requires the SeqinR package.
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
  alignment <- pairwiseAlignment(s1, s2, type = seq_align, substitutionMatrix = mat, 
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
  
  # 6. Histograma
  hist(randomscores, col="red")
  
}




seq1 <- "gi32141095_N_0.fa"
seq2 <- "gi32141095_N_1.fa"
seq_type ="dna"
seq_align = "global"
mat = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
gap = c(-2, -8)
N = 100
shuff = 1

seq_alignment(seq1, seq2, seq_type, seq_align, mat, gap, N, shuff)

library(evd)
scores <- table(randomscores)
moda <- as.numeric(names( which(scores == max(scores))[1]))
dens <- dgumbel((randomscores*(-1)), loc = moda, sd(randomscores)) 

hist(randomscores, col="grey", breaks = (N/5), freq = FALSE)
lines(dens, lwd = 3)


#### Resultado gráfico: 
  
# Obtener un histograma de los N (nº de replicas) scores obtenidos al hacer el alineamiento 
# por parejas entre la secuencia intacta versus las secuencia shuffling. 

## 1. Añadir al histograma la distribución de Gumbel que se ajusta e indicar con una flecha donde 
##    se sitúa el score original en el gráfico. 

## 2. Por otro lado, añadir dentro del gráfico como texto:
### 2.1. La información básica del alineamiento realizado (inputs): 
####     Nombre de las secuencias, 
####     tipo de alineamiento, 
####     tipo de secuencia, 
####     matriz de substitución, 
####     puntuación Gap, 
####     número de réplicas y 
####     secuencia donde se hace shuffling;  

### 2.2. Los resultados del alineamiento (outputs): 
####     parámetros de la distribución de Gumbel, 
####     score original, 
####     Probabilidad(distr Gumbel estimada > score original) y 
####     Probablidad empirica(Histograma > score original)

#### Resultado no gráfico:

## 1. De los N scores obtenidos haciendo shuffling mostrar el resultado de aplicar la función summary(). 
## 2. Por otra parte mostrar:
####   los parámetros de la distribución de Gumbel ajustada a los datos, 
####   la K, 
####   la S, 
####   la S’ y 
####   la P(S’)



