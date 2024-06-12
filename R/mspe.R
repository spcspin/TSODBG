#' The MSPE score function
#'
#' This function is designed to calculate MSPE score.
#'
#' @param X A numeric matrix of the marker-associated matrices for the training set.
#' @param X0 A numeric matrix of the marker-associated matrices for the candidate population.
#' @param lambda  An integer number of the shrinkage parameter (deafult: 1)
#'
#' @return This function will return the MSPE score value.
#'
#' @export
#' @examples
#' data(PC)
#' mspe_score(PC[1:50,],PC[51:100,])
#
mspe_score <- function(X,X0,lambda=1){
  nr <- nrow(X)
  n0r <- nrow(X0)
  n0c <- ncol(X0)
  A <- t(X)%*%solve((X%*%t(X)+diag(nr)*lambda))
  IAX <- diag(n0c) - A%*%X
  X0IAX <- X0%*%IAX
  X0A <- X0%*%A
  B <- sum(diag(X0A%*%t(X0A)))
  C <- sum(diag(t(X0IAX)%*%X0IAX))
  return(1+1/n0r*B+C)
}
