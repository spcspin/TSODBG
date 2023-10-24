#' The A-optimality function
#'
#' This function is designed to calculate A-optimality.
#'
#' @param K A numeric matrix of kinship.
#'
#' @return This function will return each individuals' A-optimality values.
#'
#' @export
#' @examples
#' data(kinship)
#' A.opt(kinship)
#
A.opt <- function(K){
  Nc <- length(K[,1])
  A.opt <- diag(K)
  ID_A.opt <- as.data.frame(cbind(A.opt,colnum=seq(1,Nc)))
  ID_A.opt$colnum <- as.numeric(ID_A.opt$colnum)
  ID_A.opt$A.opt <- as.numeric(ID_A.opt$A.opt)
  return(ID_A.opt)
}
