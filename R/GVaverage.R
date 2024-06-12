#' The GV-average function
#'
#' This function is designed to calculate GV-average.
#'
#' @param K A numeric matrix of kinship.
#'
#' @return This function will return each individuals' GV-average values.
#'
#' @export
#' @examples
#' data(kinship)
#' GVaverage(kinship)
#
GVaverage <- function(K){
  Nc <- length(K[,1])
  GVaverage <- diag(K)
  ID_GVaverage <- as.data.frame(cbind(GVaverage,colnum=seq(1,Nc)))
  ID_GVaverage$colnum <- as.numeric(ID_GVaverage$colnum)
  ID_GVaverage$GVaverage <- as.numeric(ID_GVaverage$GVaverage)
  return(ID_GVaverage)
}
