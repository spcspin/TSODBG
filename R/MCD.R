#' The MCD function
#'
#' This function is designed to calculate MCD.
#'
#' @param K A numeric matrix of kinship.
#'
#' @return This function will return each individuals' MCD values.
#'
#' @export
#' @examples
#' data(kinship)
#' MCD(kinship)
#
MCD <- function(K){
  Nc <- length(K[,1])
  MCD <- rep(NA,Nc)
  I <- diag(Nc)
  Jbar <- matrix(1/Nc,Nc,Nc)
  for (i in 1:Nc) {
    ci <- c(rep(-1/Nc,Nc))
    ci[i] <- 1-1/Nc
    MCD[i] <-( t(ci) %*% (K-1*solve((I-Jbar)+1*solve(K))) %*% ci)
  }

  ID_MCD <- as.data.frame(cbind(MCD,colnum=seq(1,Nc)))
  ID_MCD$colnum <- as.numeric(ID_MCD$colnum)
  ID_MCD$MCD <- as.numeric(ID_MCD$MCD)
  return(ID_MCD)
}
