#' The GV-overall function
#'
#' This function is designed to calculate GV-overall.
#'
#' @param K A numeric matrix of kinship.
#'
#' @return This function will return GV-overall value.
#'
#' @export
#' @examples
#' data(PC)
#' GVoverall(kinship[1:50,1:50])
#

GVoverall <- function(K){
  score <- det(K)
  return(score)
}

