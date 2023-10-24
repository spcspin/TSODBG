#' The Genetic Algorithm function
#'
#' This function is designed to perform genetic algorithm modified form Ou and Liao (2019).
#'
#' @param K A numeric matrix of kinship. (for D.opt method)
#' @param geno A numeric matrix of marker-associated matrix (for rSocre, MSPE, and D.opt methods; rows: individuals; columns: PCs or markers).
#' @param cand An integer vector of which individuals are candidates of the training set.
#' @param n.train The size of the training set.
#' @param subpop A character vector of subpopulation.
#' @param test An integer vector of which individuals are in the test set.
#' @param criteria Select the criteria uses to discover the best genotypes (deafult: r_score).
#' @param greater  A boolean variable determines whether the criteria is the greater, the better.
#' @param min.iter Minimum iteration.
#' @param console Whether the function will print out the iteration numbers (deafult: TRUE).
#'
#' @return This function will return 3 information including OPTtrain (a vector of chosen optimal training set), TOPscore (highest scores of before iteration), and ITERscore (criteria scores of each iteration).
#'
#' @export
#' @examples
#' data(PC)
#' GA(geno=PC[,1:50], cand = 1:328, n.train = 50)
#'

GA <- function(K=NULL, geno=NULL, cand, n.train, subpop=NULL, test=NULL, criteria=r_score, greater=TRUE, min.iter=NULL, console=TRUE)
{
  if (!is.null(K) && !is.null(geno)) {
    stop("You can only provide either 'K' or 'geno', not both.")
  } else if (is.null(K) && is.null(geno)) {
    stop("You need to provide either 'K' or 'geno'.")}
  n=n.train; N=nrow(geno); Nc=length(cand)
  geno = as.matrix(geno)

  if(is.null(subpop)){
    if(is.null(test)){
        ## npop untarget
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
        sol = sample(cand,n)
        if (identical(criteria,D.opt)) {
          score=criteria(K[sol, sol])
        }else{
          score = criteria(geno[sol,], geno[cand[!cand%in%sol],])
        }

        if (greater==FALSE) {
          score=-score
        }
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          new.sol[sample(n,1)] = sample(cand[!cand%in%sol],1)
          if (identical(criteria,D.opt)) {
            new.score=criteria(K[new.sol, new.sol])
          }else{
            new.score=criteria(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          }
          if (greater==FALSE) {
            new.score=-new.score
          }
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score, new.score)
            top.score = c(top.score, new.score)
          }else{
            iter.score=c(iter.score, new.score)
            top.score=c(top.score, score)
          }
          if(console) cat(iter, "..", sep="")
          iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }

    }else{
        ## npop target
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
        sol=matrix(NA, n, 30); for(i in 1:30){sol[,i] = sample(cand, n)}
        if (identical(criteria,D.opt)) {
          score=rep(NA,30); for(i in 1:30){score[i]=criteria(K[sol[,i], sol[,i]])}
        }else{
          score=rep(NA,30); for(i in 1:30){score[i]=criteria(geno[sol[,i],], geno[test,])}
        }
        if (greater==FALSE) {
          score=-score
        }
        iter.score=score; top.score=min(score)
        stop=0; iter=1
        while(stop==0){
          elite = which(rank(score)<4)
          del = sample(seq(30)[-elite], 3, prob=(score[-elite]/sum(score[-elite])))
          for(i in 1:length(del)){
            par = sample(seq(30)[-del],2)
            sol[,del[i]] = sample(unique(c(sol[,par[1]],sol[,par[2]])),n)
          }
          for(i in 1:30){
            new.sol=sol[,i]
            if(i %in% del){
              sol[sample(n,ceiling(n*.03)),i]=sample(cand[!cand%in%sol[,i]],ceiling(n*.03))
              if (identical(criteria,D.opt)) {
                score[i] = criteria(K[sol[,i], sol[,i]])
              }else{
                score[i] = criteria(geno[sol[,i],], geno[test,])
              }
              if (greater==FALSE) {
                score[i]=-score[i]
              }
            }else{
              new.sol[sample(n,ceiling(n*.03))]=sample(cand[!cand%in%new.sol],ceiling(n*.03))
              if (identical(criteria,D.opt)) {
                score[i] = criteria(K[new.sol,new.sol])
              }else{
                new.score = criteria(geno[new.sol,], geno[test,])
              }

              if (greater==FALSE) {
                new.score=-new.score
              }
              if(new.score > score[i]){sol[,i] = new.sol; score[i]=new.score}
            }
          }
          iter.score=c(iter.score,mean(score)); top.score=c(top.score, min(score))
          if(console) cat(iter, "..", sep="")
          iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }
        sol = sol[,which(score==min(score))[1]]

    }
  }else{
    subpop=as.character(subpop)
    if(is.null(test)){
        ## pop untarget
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
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
          sol = c(sol, sample(cand[subpop==pops[i]],pop.ratio[i]))
        }
        if (identical(criteria,D.opt)) {
          score=criteria(K[sol, sol])
        }else{
          score = criteria(geno[sol,], geno[cand[!cand%in%sol],])
        }
        if (greater==FALSE) {
          score=-score
        }
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1)
          new.sol[sample(which(subpop[sol]==pops[p]),1)]=sample(which(subpop==pops[p])[!which(subpop==pops[p])%in%sol],1)
          if (identical(criteria,D.opt)) {
            new.score=criteria(K[new.sol, new.sol])
          }else{
            new.score = criteria(geno[new.sol,], geno[cand[!cand%in%new.sol],])
          }
          if (greater==FALSE) {
            new.score=-new.score
          }
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,new.score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          if(console) cat(iter, "..", sep="")
          iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }


    }else{
        ## pop target pev
        if(is.null(min.iter)){min.iter=round(sqrt(Nc*n)*50)}
        pops = names(table(subpop))
        pop.ratio = ceiling(n*as.numeric(table(subpop))/sum(as.numeric(table(subpop))))
        stop=0
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else{
            stop=1
          }
        }
        sol=c()
        for(i in 1:length(pops)){
          sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
        }
        if (identical(criteria,D.opt)) {
          score=criteria(K[sol, sol])
        }else{
          score=criteria(geno[sol,], geno[test,])
        }
        if (greater==FALSE) {
          score=-score
        }
        iter.score=score; top.score=score
        stop=0; iter=1
        while(stop==0){
          new.sol=sol
          p = sample(length(pops),1,prob=pop.ratio)
          new.sol[as.numeric(sample(as.character(which(subpop[new.sol]==pops[p])),1))]=cand[sample(which(subpop==pops[p])[!which(subpop==pops[p])%in%sol],1)]
          if (identical(criteria,D.opt)) {
            new.score=criteria(K[new.sol, new.sol])
          }else{
            new.score=criteria(geno[new.sol,],geno[test,])
          }

          if (greater==FALSE) {
            new.score=-new.score
          }
          if(new.score > score){
            sol=new.sol; score=new.score
            iter.score=c(iter.score,score)
            top.score=c(top.score,score)
          }else{
            iter.score=c(iter.score,new.score)
            top.score=c(top.score,score)
          }
          if(console) cat(iter, "..", sep="")
          iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-sqrt(Nc*n)*5])<1e-6){stop=1}}
          if(iter > (min.iter*2)){stop=1}
        }

    }
  }
  ret = list(
    OPTtrain=as.numeric(sol),
    TOPscore=as.numeric(top.score[-1]),
    ITERscore=as.numeric(iter.score[-1])
  )
  return(ret)
}
