#' INTERNAL: Converts a matrix with dummies in columns into categorical vector
#'
#' @param D a matrix of dummies.
#' @author Maciej J. Danko
#' @keywords internal
DummyMat2Vector<-function(D) D %*% ((1L : dim(D)[2L]) -1L)


contingencytables <- function(model, formula, data, names.reg=identity){
  if (class(names.reg)=='function') NN <- names.reg(colnames(model$reg.mm)) else NN <- names.reg
  if (class(formula)=='formula')
    tmp <- formula2classes(formula, data, sep=' ', return.matrix = TRUE) else stop('The formula parameter must be of class "formula".')
  colnames(model$reg.mm) <- NN
  M <- tmp$x
  cTAB <- sapply(levels(M), function(k) colSums(model$reg.mm[k==M,]))
  cTAB <- cbind(cTAB, ALL=rowSums(cTAB))
  fTAB.1 <- 100*cTAB / length(M)
  cTAB <- rbind(cTAB, 'Number of subjects'=c(table(tmp$x),length(M)))
  #fTAB.2 <- (100 * cTAB / t(matrix(cTAB[dim(cTAB)[1],],dim(cTAB)[2],dim(cTAB)[1])))[-dim(cTAB)[1],]
  fTAB.2 <- (100 * cTAB / matrix(cTAB[,dim(cTAB)[2]],dim(cTAB)[1],dim(cTAB)[2]))[,-dim(cTAB)[2]]
  list(countsTAB=cTAB, freqTAB1=fTAB.1, freqTAB2=fTAB.2)
}
