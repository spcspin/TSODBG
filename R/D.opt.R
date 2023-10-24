#' The D-optimality function
#'
#' This function is designed to calculate D-optimality.
#'
#' @param K A numeric matrix of kinship.
#'
#' @return This function will return D-optimality value.
#'
#' @export
#' @examples
#' data(PC)
#' D.opt(kinship[1:50,1:50])
#

D.opt <- function(K){
  score = det(K)
  return(score)
}
