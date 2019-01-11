#' @title Quantile regression using R.
#' @description Quantile regression by Linear Programming using R. Which means turning Quantile regression into  a linear programming problem
#' @import lpSolve
#' @param x the observed value of independent variables x should be a vector or matrix.(numeric)
#' @param y the observed value of dependent variable and y should be a vector. (numeric)
#' @param tau the quantile level and tau should be in (0,1). (numeric)
#' @param Intercept if TRUE then the reslut including the intercept coefficient (logical)
#' @return the regression coefficients
#' @examples
#' \dontrun{
#' data(engel)
#' b1hat<-matrix(NA,ncol=2,nrow=6);b2hat<-matrix(NA,ncol=2,nrow=6)
#' taus <- c(.05,.1,.25,.75,.9,.95)
#' for(i in 1:length(taus)){
#'  b1hat[i,]=qrlinear(income,foodexp,taus[i],Intercept=T)
#'  b2hat[i,]=rq(foodexp~income,taus[i])$coef
#'  }
#'  print(b1hat);print(b2hat)
#' }
#' @export
qrlinear<-function(x,y,tau,Intercept=T){
  if(tau>=1 | tau<=0) stop("tau is not in (0,1)")
  else{
    if(Intercept) X=cbind(1,x)
    else          X=as.matrix(x)
    n=dim(X)[1];p=dim(X)[2]
    obj <- c(rep(0,2*p),rep(1-tau,n),rep(tau,n))
    dir <- rep("<=",2*n); rhs <- c(y,-y)
    con<-cbind(rbind(X,-X),rbind(-X,X),-diag(2*n,2*n))
    res.lp <- lp("min",obj,con,dir,rhs)
    b1hat <- res.lp$solution[1:p]-res.lp$solution[p+1:p]
  }
  if(Intercept)  {
    names=paste(c("x"), 1:(p-1),sep="")
    names=c("Intercept",names)
    names(b1hat)=names
  }
  else{
    names=paste(c("x"), 1:p,sep="")
    names(b1hat)=names
  }
  return(b1hat)
}
