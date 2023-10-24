#' The rScore function
#'
#' This function is designed to calculate rScore from Ou and Liao (2019).
#'
#' @param X A numeric matrix of the marker-associated matrices for the training set.
#' @param X0 A numeric matrix of the marker-associated matrices for the candidate population.
#'
#' @return This function will return the rScore value.
#'
#' @export
#' @examples
#' data(PC)
#' r_score(PC[1:50,],PC[51:100,])
#
r_score <- function(X,X0){
  nr = nrow(X)
  nc = ncol(X)
  n0r = nrow(X0)
  A = t(X)%*%solve((X%*%t(X)+diag(nr)*1/nc))
  IJ = diag(n0r)-(diag(n0r)*1/n0r)
  q1 = n0r-1 + sum((IJ %*% X0)^2)
  IJX0A = IJ%*%X0%*%A
  IJX0AX = IJX0A%*%X
  q2 = sum(IJX0A^2) + sum(IJX0AX^2)
  q12 = sum(diag(t(X0)%*%IJX0AX))
  return(q12/sqrt(q1*q2))
}
