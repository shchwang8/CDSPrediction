codonList <- c(
  "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG",
  "AAT", "AAC", "GAT", "GAC", "TGT", "TGC", "CAA", "CAG", "GAA", "GAG",
  "GGT", "GGC", "GGA", "GGG", "CAT", "CAC", "ATT", "ATC", "ATA", "TTA",
  "TTG", "CTT", "CTC", "CTA", "CTG", "AAA", "AAG", "ATG", "TTT", "TTC",
  "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
  "ACT", "ACC", "ACA", "ACG", "TGG", "TAT", "TAC", "GTT", "GTC", "GTA",
  "GTG", "TAA", "TGA", "TAG")

  # String Find: pos <- instr(fnm, "_") + 1
  instr <- function(str1, str2, startpos=1, n=1){
    aa=unlist(strsplit(substring(str1, startpos), str2))
    if(length(aa) < n+1 ) return(0);
    return(sum(nchar(aa[1:n])) + startpos + (n-1)*nchar(str2) )
  }#instr

  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  }#mkdirs

  filestem <- function(filenm){
    len <- nchar(filenm)
    str1 <- substr(filenm, len-3, len)
    if(str1 == ".txt") return(substr(filenm, 1, len-4))
  }#filestem

#' Create a file with Codon count numbers
#'
#' @param infile DNA sequence file
#' @return None (new file with "*_codonSEQ.txt" will be created)
#' @export
  CodonCount <- function(filenm) {
    fnm <- filenm
    DD <- read.table(fnm, header = TRUE, sep = "\t")
    N = nrow(DD)
    M <- length(codonList)
    colnames(DD)[1] <- "Gene_ID"
    DD$RanVal <- runif(N)*10
    for(c in c(1:M)){
      DD[c + 3] <- rep(0, N)
      colnames(DD)[c + 3] <- codonList[c]
    }
    DD <- DD[, colnames(DD)[c(1, 3, c(4:(M+3)), 2)]]
    M <- ncol(DD)
    i <- 1
    for(i in c(1:N)){
      x <- toupper(toString(DD[i, M]))
      strVec <- substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3))
      ctTab <- as.data.frame(table(strVec))
      for(j in c(3:(M-1))){
        nm <- colnames(DD)[j]
        val <- ctTab[ctTab$strVec==nm, 2]
        if(length(val) > 0) DD[i, j] <- val
      } #for j
      DD[i,]
    } #for i
    fnm0 <- filestem(filenm)
    fnm <- paste(fnm0, "_codonSEQ.txt", sep = "")
    write.table(DD, file = fnm, quote = FALSE, sep = "\t", row.names = FALSE)
  }#CodonCount

#' Use build models to predict CDS expression value
#'
#' @param infile DNA sequence file
#' @return None
#' @export
  autoMLpredict<- function(filenm){
    dFn <- paste(filestem(filenm), "_codonSEQ.txt", sep = "")
    dD <- h2o.importFile(path = normalizePath(dFn))
    # dim(dD)
    m <- 1
    for(m in c(1:23)){
      mFn <- paste(tPath, "/", modelFD, "/T", m, "_AutoML_2020", sep = "")
      m2 <- h2o.loadModel(mFn)
      # m2
      predictions <- h2o.predict(object = m2, newdata = dD)
      pred <- as.data.frame(predictions$predict)
      dDD <- as.data.frame(dD)
      dDD$pred <- as.vector(pred[[1]])
      M <- ncol(dDD)
      dDD <- dDD[, colnames(dDD)[c(1, M, c(2:(M-1)))]]
      fnm <- paste(resultFD, "/", filenm, "_T", m, "-coeResult.txt", sep = "")
      write.table(dDD, fnm, row.names = FALSE, sep = "\t")
    } #m
  }#autoMLpredict

#' Get the summary data file
#'
#' @param infile DNA sequence file
#' @return None
#' @export
  summaryPredict <- function(filenm){
    m <- 1
    fnm <- paste(resultFD, "/", filenm, "_T", m, "-coeResult.txt", sep = "")
    DD <- read.table(fnm, header = TRUE, sep = "\t")
    M <- ncol(DD)
    DD <- DD[, colnames(DD)[c(1, c(4:M))]]
    MM <-23
    for(m in c(1:MM)){
      fnm <- paste(resultFD, "/", filenm, "_T", m, "-coeResult.txt", sep = "")
      dD <- read.table(fnm, header = TRUE, sep = "\t")
      DD$pred <- dD[, 2]
      M <- ncol(DD)
      clNm <- paste("AutoML_T", m, sep = "")
      colnames(DD)[M] <- clNm
    } #m
    colnames(DD)
    M <- ncol(DD)
    DD <- DD[, colnames(DD)[c(1, c((M-MM+1):M), c(2:(M-MM)))]]
    colnames(DD)
    numV <- DD[, 2:(MM+1)]
    mean <- rowMeans(numV)
    DD$AveAllTissues <- mean
    colnames(DD)
    M <- ncol(DD)
    DD <- DD[, colnames(DD)[c(1, M, c(2:(M-1)))]]
    colnames(DD)
    numV <- DD[, 14:17]
    mean <- rowMeans(numV)
    DD$AveLeaf <- mean
    colnames(DD)
    M <- ncol(DD)
    DD <- DD[, colnames(DD)[c(1, M, c(2:(M-1)))]]
    colnames(DD)
    newID <- paste(filestem(filenm), "_", DD[, 1],sep = "")
    DD$ID <- newID
    colnames(DD)
    M <- ncol(DD)
    DD <- DD[, colnames(DD)[c(M, c(1:(M-1)))]]
    colnames(DD)
    fonm <- paste(sumFD, "/R", filestem(filenm), "-Summary.txt", sep = "")
    write.table(DD, fonm, quote = FALSE, row.names = FALSE, sep = "\t")
  }#summaryPredict


