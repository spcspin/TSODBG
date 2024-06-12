#' The CDmean(v2) function
#'
#' This function is designed to calculate CDmean(v2).
#'
#' @param K A numeric matrix of kinship.
#' @param train An integer vector of which individuals are in the training set.
#' @param test An integer vector of which individuals are in the test set.
#'
#' @return This function will return the CDmean(v2) value.
#'
#' @export
#' @examples
#' data(kinship)
#' CDmeanv2(kinship,1:50,1:100)
#
CDmeanv2 <- function(K, train, test){
  Kt <- K[train,train]
  K0 <- K[test,test]
  Kt0 <- K[train,test]
  K0t <- K[test,train]
  n0 <- length(test)
  nt <- length(train)
  I <- diag(nt)
  Jbar <- matrix(1/nt,nt,nt)
  M <- I-Jbar
  Bi <- diag(K0t%*%solve(M%*%Kt+1*I)%*%M%*%t(K0t))
  Ai <- diag(K0)
  CDv2 <- sum(Bi/Ai)
  return(CDv2)
}
