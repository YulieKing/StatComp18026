#' @title Linear regression using R.
#' @description Linear regression  by gradient descent using R.
#' @importFrom stats rnorm
#' @param x The observed value of independent variables x should be a vector or matrix.(numeric)
#' @param y The observed value of dependent variable and y should be a vector. (numeric)
#' @param crit The value that controls convergence. (numeric)
#' @param maxiter Maximum number of iteration steps. (numeric)
#' @param Intercept if TRUE then the reslut including the intercept coefficient (logical)
#' @return the regression coefficients
#' @examples
#' \dontrun{
#' x1=rnorm(100);x2=rnorm(100)
#' x=cbind(x1,x2)
#' theta=c(1.1,2)
#' y=x%*%theta+3+rnorm(100)
#' GDLM(x,y,crit=1e-12,500,Intercept=T)
#' }
#' @export
GDLM<-function(x,y,crit=1e-8,maxiter,Intercept=T){
  iter=1;newcrit=1
  if(Intercept) x=cbind(1,as.matrix(x))
  else x=as.matrix(x)
  m=nrow(x);n=ncol(x);theta=matrix(rep(0,n),n,1)
  while((newcrit>crit) | (iter<maxiter)){
    iter=iter+1
    h=x%*%theta
    des=t(t(h-y)%*%x)
    sstep=1;alpha=0.25;beta=0.7
    new_theta=theta-sstep*des
    new_h=x%*%new_theta
    costfunction=t(h-y)%*%(h-y)
    new_costfunction=t(new_h-y)%*%(new_h-y)
    while(new_costfunction>costfunction-alpha*sstep*sum(des*des)){
      sstep=sstep*beta
      new_theta=theta-sstep*des
      new_h=x%*%new_theta
      new_costfunction=t(new_h-y)%*%(new_h-y)
    }
    newcrit<-t(theta-new_theta)%*%(theta-new_theta)
    theta=new_theta
  }
  if(Intercept)  {
    names=paste(c("x"), 1:(n-1),sep="")
    names=c("Intercept",names)
    rownames(theta)=names
  }
  else{
    names=paste(c("x"), 1:n,sep="")
    rownames(theta)=names
  }
  costfunction=t(x%*%theta-y)%*%(x%*%theta-y)
  result=list(theta,iter,costfunction)
  names(result)<-c('coef','iter','error')
  return(result)
}
