#' @title generate random variables from a Cauchy distribution using R.
#' @description Use the Metropolis-Hastings sampler to generate random variables from a Cauchy distribution.
#' @importFrom stats dcauchy
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @param n length of the Markov chain(numeric)
#' @param sigma standard error of the proposal distribution (numeric)
#' @param x0 initial value (numeric)
#' @param u parameter of cauchy distribution. (numeric)
#' @param lamada parameter of cauchy distribution. (numeric)
#' @return the Markov chain and the number of refuse.
#' @examples
#' \dontrun{
#' n<-10000
#' cauchy<-rdCauchy(n,sigma=0.5,x0=0,u=0,lamada=1)
#' refuse.pr<-cauchy$k/n
#' refuse.pr
#' hist(cauchy$x[1001:n],freq=F,main="density of cauchy",breaks=60)
#' curve(dcauchy(x,0,1),add=TRUE)
#' }
#' @export
rdCauchy<-function(n,sigma,x0,u,lamada){
  x<-NULL
  x[1]<-x0
  e=runif(n)
  k<-0
  for(i in 2:n){
    y<-rnorm(1,x[i-1],sigma)
    if(e[i]<=(dcauchy(y,u,lamada)/dcauchy(x[i-1],u,lamada)))  x[i]<-y
    else{
      x[i]<-x[i-1]
      k<-k+1
    }
  }
  return(list(x=x,k=k))
}

