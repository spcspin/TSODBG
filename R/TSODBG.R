#' Training Set Optimization to Discover the Best Genotypes
#'
#' This function is designed for discover the best genotypes.
#'
#' @param K A numeric matrix of kinship. (for D.opt, A.opt, CD, and MCD methods)
#' @param geno A numeric matrix of marker-associated matrix (for rSocre, MSPE, and D.opt methods; rows: individuals; columns: PCs or markers).
#' @param cand An integer vector of which individuals are candidates of the training set.
#' @param n.train The size of the training set.
#' @param subpop A character vector of subpopulation.
#' @param test An integer vector of which individuals are in the test set.
#' @param method Select the method uses to discover the best genotypes (deafult: rScore).
#' @param min.iter Minimum iteration.
#' @param console Whether the function will print out the iteration numbers (deafult: TRUE).
#'
#' @return This function will return the best genotypes.
#'
#' @export
#' @examples
#' data(PC)
#' TSODBG(geno=PC[,1:50], cand = 1:328, n.train = 50)
#'

TSODBG<-function(K=NULL, geno=NULL, cand, n.train, subpop=NULL, test=NULL, method="rScore", min.iter=NULL, console=TRUE)
{
  if (!is.null(K) && !is.null(geno)) {
    stop("You can only provide either 'K' or 'geno', not both.")
  } else if (is.null(K) && is.null(geno)) {
    stop("You need to provide either 'K' or 'geno'.")
  } else if (!is.null(K)) {
    if (method=="A.opt"){
      n=n.train; ID_A.opt=A.opt(K)
      if (is.null(subpop)) {
        sol=sort(ID_A.opt[order(ID_A.opt[,1],decreasing = T)[1:n.train],2])
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
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, order(ID_A.opt[,1],decreasing = T)[subpop==pops[i]][1:pop.ratio[i]] )
          sol = sort(sol)
        }
      }
    }
    if (method=="CD"){
      n=n.train; ID_CD=CD(K)
      if (is.null(subpop)) {
        sol=sort(ID_CD[order(ID_CD[,1],decreasing = T)[1:n.train],2])
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
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, order(ID_CD[,1],decreasing = T)[subpop==pops[i]][1:pop.ratio[i]] )
          sol = sort(sol)
        }
      }
    }
    if (method=="MCD"){
      n=n.train; ID_MCD=MCD(K)
      if (is.null(subpop)) {
        sol=sort(ID_MCD[order(ID_MCD[,1],decreasing = T)[1:n.train],2])
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
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, order(ID_MCD[,1],decreasing = T)[subpop==pops[i]][1:pop.ratio[i]] )
          sol = sort(sol)
        }
      }
    }
  } else {
    if (method=="rScore") {
      sol=GA(K=NULL,geno,cand,n.train, subpop, test, greater = T, criteria = r_score, min.iter, console)$OPTtrain
      sol = sort(sol)
    }
    if (method=="MSPE") {
      sol=GA(K=NULL,geno,cand,n.train, subpop, test, greater = F, criteria = mspe_score, min.iter, console)$OPTtrain
      sol = sort(sol)
    }
    if (method=="D,opt") {
      sol=GA(K, geno=NULL,cand,n.train, subpop, test, greater = T, criteria = D.opt, min.iter, console)$OPTtrain
      sol = sort(sol)
    }
  }
  return(sol)

}

