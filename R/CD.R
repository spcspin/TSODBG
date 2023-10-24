#' The CD function
#'
#' This function is designed to calculate CD.
#'
#' @param K A numeric matrix of kinship.
#'
#' @return This function will return each individuals' CD values.
#'
#' @export
#' @examples
#' data(kinship)
#' CD(kinship)
#
CD <- function(K){
  Nc <- length(K[,1])
  CD <- rep(NA,Nc)
  I <- diag(Nc)
  Jbar <- matrix(1/Nc,Nc,Nc)
  for (i in 1:Nc) {
    ci <- c(rep(-1/Nc,Nc))
    ci[i] <- 1-1/Nc
    CD[i] <-( t(ci) %*% (K-1*solve((I-Jbar)+1*solve(K))) %*% ci) / (t(ci) %*% K %*% ci)
  }

  ID_CD <- as.data.frame(cbind(CD,colnum=seq(1,Nc)))
  ID_CD$colnum <- as.numeric(ID_CD$colnum)
  ID_CD$CD <- as.numeric(ID_CD$CD)
  return(ID_CD)
}
