#' Training Set Optimization to Discover the Best Genotypes
#'
#' This function is designed for discover the best genotypes.
#'
#' @param K A numeric matrix of kinship. (for  GVoverall, GVaverage, CDmean(v2), and CDranking methods)
#' @param geno A numeric matrix of marker-associated matrix or principal components (for rSocre and MSPE methods; rows: individuals; columns: PCs or markers).
#' @param cand An integer vector of which individuals are candidates of the training set.
#' @param n.train The size of the training set.
#' @param subpop A character vector of subpopulation.
#' @param test An integer vector of which individuals are in the test set.
#' @param method Select the method uses to discover the best genotypes (deafult: rScore).
#' @param min.iter Minimum iteration.
#' @param max.iter Maximum iteration.
#' @param console Whether the function will print out the iteration numbers (deafult: TRUE).
#'
#' @return This function will return the best genotypes.
#'
#' @export
#' @examples
#' data(PC)
#' TSODBG(geno=PC[,1:50], cand = 1:328, n.train = 50, min.iter = 50)
#'

TSODBG<-function(K=NULL, geno=NULL, cand, n.train, subpop=NULL, test=NULL, method="rScore", min.iter=NULL, max.iter=NULL, console=TRUE)
{
  if (!is.null(K) && !is.null(geno)) {
    stop("You can only provide either 'K' or 'geno', not both.")
  } else if (is.null(K) && is.null(geno)) {
    stop("You need to provide either 'K' or 'geno'.")
  } else if (!is.null(K)) {
    if (method=="GVaverage"){
      n=n.train; ID_GVaverage=GVaverage(K)
      if (is.null(subpop)) {
        result.sol=sort(ID_GVaverage[order(ID_GVaverage[,1],decreasing = T)[1:n.train],2])
      }else{
        if(length(cand)!=length(subpop)){stop("Input data is not correct.")}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*(as.numeric(table(subpop))/sum(as.numeric(table(subpop)))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{stop=1}
        }
        result.sol=c()
        for(i in 1:length(pops)){
          result.sol = c(result.sol, order(ID_GVaverage[,1],decreasing = T)[subpop==pops[i]][1:pop.ratio[i]] )
          result.sol = sort(result.sol)
        }
      }
    }
    if (method=="CDranking"){
      n=n.train; ID_CDranking=CDranking(K)
      if (is.null(subpop)) {
        result.sol=sort(ID_CDranking[order(ID_CDranking[,1],decreasing = T)[1:n.train],2])
      }else{
        if(length(cand)!=length(subpop)){stop("Input data is not correct.")}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*(as.numeric(table(subpop))/sum(as.numeric(table(subpop)))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{stop=1}
        }
        result.sol=c()
        for(i in 1:length(pops)){
          result.sol = c(result.sol, order(ID_CDranking[,1],decreasing = T)[subpop==pops[i]][1:pop.ratio[i]] )
          result.sol = sort(result.sol)
        }
      }
    }

  } else {
    if (method=="rScore") {
      if(is.null(max.iter)){max.iter=10000}
      result.sol = GA(K=NULL,geno,cand,n.train, subpop, test, greater = T, criteria = r_score, min.iter, max.iter, console)$OPTtrain
      result.sol = sort(result.sol)
    }
    if (method=="MSPE") {
      if(is.null(max.iter)){max.iter=10000}
      result.sol = GA(K=NULL,geno,cand,n.train, subpop, test, greater = F, criteria = mspe_score, min.iter, max.iter, console)$OPTtrain
      result.sol = sort(result.sol)
    }
    if (method=="GVoverall") {
      if(is.null(max.iter)){max.iter=10000}
      result.sol = GA(K, geno=NULL,cand,n.train, subpop, test, greater = T, criteria = GVoverall, min.iter, max.iter, console)$OPTtrain
      result.sol = sort(result.sol)
    }
    if (method=="CDmeanv2") {
      if(is.null(max.iter)){max.iter=25000}
      result.sol = GA(K, geno=NULL,cand,n.train, subpop, test, greater = T, criteria = CDmeanv2, min.iter, max.iter, console)$OPTtrain
      result.sol = sort(result.sol)
    }
  }
  return(result.sol)

}

