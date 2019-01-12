## ------------------------------------------------------------------------
A<-matrix(1:9,3,3)
B<-diag(2,3,3)
C<-A%*%B
C

## ------------------------------------------------------------------------
library(lattice)
n <- seq(5, 20, 5)
x <- rnorm(sum(n))
y <- factor(rep(n, n), labels=paste("n =", n))
densityplot(~ x | y,
panel = function(x, ...) {
panel.densityplot(x, col="DarkOliveGreen", ...)
panel.mathdensity(dmath=dnorm,
args=list(mean=mean(x), sd=sd(x)),
col="darkblue")
})

## ------------------------------------------------------------------------
ts(1:47, frequency = 12, start = c(2018, 9))

## ------------------------------------------------------------------------
set.seed(123)
#The inverse transform method
n=1000                        #size of sample
u=runif(n)                    #generate u from [0,1] uniform distribution
x=0:4                         #The set of values of variable X 
p=c(.1,.2,.2,.2,.3)           #The probability of X
cp=cumsum(p)                  #The cumulative probability of X
r=x[findInterval(u,cp)+1]     #Find the corresponding value of x
ct <- as.vector(table(r))     #Compare to the real probability
ct/sum(ct)/p                  #A value close to 1 indicates a more accurate result

#Repeat using the R sample function.
r.sample<- sample(x, size = 1000, replace = TRUE,prob = p)
ct.sample <- as.vector(table(r.sample))
ct.sample/sum(ct.sample)/p     

## ------------------------------------------------------------------------
random.Beta<-function(a,b,n){
 if(a<=0||b<=0) stop("alphas and betas must be positive")
 else if(a<=1||b<=1)
   stop("The maximum of Beta(a,b) is infinite,so This method is not applicable")
 else {
   j<-k<-0;y <- numeric(n)
   while (k < n) {
   u <- runif(1)
   j <- j + 1
   x <- runif(1)      #random variate from g
   x0=(a-1)/(a+b-2)   #x0 is the maximum point of the theoretical density of Beta(a,b)
   c=dbeta(x0,a,b)    #c is the maximum of the theoretical density of Beta(a,b)
   if (x * (1-x) > c*u) {             #we accept x
     k <- k + 1
     y[k] <- x
     }
   }
}
 return(y)
}
x=random.Beta(3,2,1000)
hist(x, prob = T, main ="Histogram of random observations
from Beta(3,2)   vs.  The theoretical density of Beta(3,2)",breaks=seq(0,1,0.05))
curve(dbeta(x,2,2),add=T,col="red")
#random.Beta(3,2,1000)
    ##jie suan fa 

## ------------------------------------------------------------------------

pdf.Y<-function(x,r,beta){
      y=r*beta^r/((beta+y)^(r+1))
      return(y)
}
n <- 1000 
r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda)
hist(x, prob = T, main ="Histogram of random observations
from this mixture    vs.  The theoretical density of Y",breaks=seq(0,12,0.5))
y <- seq(0, 12, .1)
lines(y,pdf.Y(y,r,beta),col="red")  ## red line is the theoretical density of Y

## ------------------------------------------------------------------------
set.seed(15656)
MC.Beta33<-function(x){
   m=10000     # The number of random Numbers
   n=length(x) 
   y=NULL
   if(x<=0||x>=1) stop("x must be in (0,1)")
   else {
     for (i in 1:n) {
        g=x[i]*dbeta(runif(m,min=0,max=x[i]),3,3)  #Generate MC random Numbers
        y[i]=mean(g) #the  estimated values
     }
   }
 return(y)
}
x=seq(.1,.9,.1)
ture=pbeta(x,3,3) 
y=MC.Beta33(x) 
print(ture)   ##print the values returned by the pbeta function in R
print(y)      ##print the value of simulation by MC method
print(y/ture) ##print the ratio

## ------------------------------------------------------------------------
library(VGAM)
MC.Rayleigh<-function(x,sigma,antithetic = TRUE){
   m=10000     # The number of random Numbers
   n=length(x) 
   y=NULL
   if(x<0||sigma<=0) stop("x should be non-negative and sigma  should be positive")
   else {
     for (i in 1:n) {
        #Using antithetic variables to generate MC random Numbers
        u=runif(m/2,min=0,max=x[i])
        if (!antithetic) v <- runif(m/2,min=0,max=x[i]) else v <- x[i] - u
        U=c(u,v)
        g=x[i]*drayleigh(U, sigma)
        y[i]=mean(g)  #the  estimated values
     }
   }
 return(y)
}
##using the rayleigh distribution function in R to check the correctness of the estimation result 
y1=MC.Rayleigh(0.3,1,antithetic = F)  # the  estimated values without using antithetic variables
y2=MC.Rayleigh(0.3,1)                 # the  estimated values by using antithetic variables
ture=prayleigh(0.3, 1)
print(c(y1/ture,y2/ture))             # The ratio shows the outcome is ture.

## Computing the percent reduction 
m=1000
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
  MC1[i] <- MC.Rayleigh(1, 2, anti = FALSE)
  MC2[i] <- MC.Rayleigh(1, 2)
}
c(sd(MC1),sd(MC2),sd(MC2)/sd(MC1))

## ------------------------------------------------------------------------
x <- seq(1, 3, .1)
w <- 2
g <- x^2/sqrt(2*pi)*exp(-x^2/2)
f1 <- exp(-x^2/2)/sqrt(2*pi)
f2<- x^2*exp(-x/2)/16
gs <- c(expression(g(x)==x^2/sqrt(2*pi)*e^{-x^2/2}),expression(f[1](x)==e^{-x^2/2}/sqrt(2*pi)),expression(f[2](x)==e^{-x/2}*{x^2}/ 16))
 par(mfrow=c(1,2))
    #figure (a)
    plot(x, g, type = "l", ylab = "",
         ylim = c(0,1), lwd = w,col=1,main='(A)')
    lines(x,f1, lty = 2, lwd = w,col=2)
    lines(x,f2, lty = 3, lwd = w,col=3)
    legend("topleft", legend = gs,
           lty = 1:3, lwd = w, inset = 0.02,col=1:3,cex=0.7)

    #figure (b)
    plot(x, g/f1, type = "l", ylab = "",
        ylim = c(0,10), lwd = w, lty = 2,col=2,main='(B)')
    lines(x, g/f2, lty = 3, lwd = w,col=3)
    legend("topleft", legend = gs[-1],
           lty = 2:3, lwd = w, inset = 0.02,col=2:3,cex=0.7)

## ------------------------------------------------------------------------
g.cut<-function(x,t){   # t is the lower limit of integral
   n=length(x)
   y=numeric(n)
   for(i in 1:n) {
     if(x[i]>t) y[i]=x[i]^2/sqrt(2*pi)*exp(-x[i]^2/2)
     else y[i]=0
   }
   return(y)
}
MC.g<-function(t){ # t is the lower limit of integral
   m=1000000     # The number of random Numbers
   n=length(t) 
   y=numeric(n)
     for (i in 1:n) {
        u=rchisq(m,6) 
        G <- g.cut(u,t[i])
        f2 <- u^2*exp(-u/2)/16
        y[i]=mean(G/f2) #the  estimated values
     }
  return(y)
} 
MC.g(0) ## this outcome shows my program is ture.
MC.g(1) ## the  estimated value

## ------------------------------------------------------------------------
set.seed(15411)
library(MASS)
Gini<-function(x){
  n<-length(x)
  a<-seq(1-n,n-1,2)
  x.i<-sort(x)
  Gini.hat<-sum(a*x.i)/(n*n*mean(x))
  return(Gini.hat)
} # you can estimate a G.hat if there comes a sample


## X is standard lognormal
n=200
m<-2000
Gini.hat1<-numeric(m)
for(i in 1:m){
  y<-rnorm(n)
  x<-exp(y) # then x is standard lognormal
  Gini.hat1[i]<-Gini(x)
}
result1<-c(mean(Gini.hat1),quantile(Gini.hat1,probs=c(0.5,0.1)))
names(result1)<-c("mean","median","deciles")
print(result1)
hist(Gini.hat1,breaks=seq(min(Gini.hat1)-0.01,max(Gini.hat1)+0.01,0.01),freq =F,main = "Histogram of Gini",xlab = "Standard lognormal")

## X is uniform
Gini.hat2<-numeric(m)
for(i in 1:m){
  x<-runif(n) # then x is uniform
  Gini.hat2[i]<-Gini(x)
}
result2<-c(mean(Gini.hat2),quantile(Gini.hat2,probs=c(0.5,0.1)))
names(result2)<-c("mean","median","deciles")
print(result2)
hist(Gini.hat2,breaks =seq(min(Gini.hat2)-0.01,max(Gini.hat2)+0.01,0.01) ,freq =F,main = "Histogram of Gini",xlab = "Uniform")

## x is  Bernoulli(0.1)
Gini.hat3<-numeric(m)
for(i in 1:m){
  x<-rbinom(n,1,0.1) # then x is Bernoulli(0.1)
  Gini.hat3[i]<-Gini(x)
}
result3<-c(mean(Gini.hat3),quantile(Gini.hat3,probs=c(0.5,0.1)))
names(result3)<-c("mean","median","deciles")
print(result3)
hist(Gini.hat3,breaks=seq(min(Gini.hat3)-0.01,max(Gini.hat3)+0.01,0.01),freq =F,main = "Histogram of Gini",xlab = "Bernoulli(0.1)")


## ------------------------------------------------------------------------
n<-200
m<-2000
mu<-0
sigma<-10
Gini1<-function(n,mu,sigma){
  y<-rnorm(n,mu,sigma)
  x<-exp(y)
  Gini.sample<-Gini(x)
  return(Gini.sample)
}
Gini.sp<-numeric(m)
for(i in 1:m){
Gini.sp[i]<-Gini1(n,mu,sigma)
}
CI<-c(mean(Gini.sp)-sd(Gini.sp)*qt(0.975,m-1),mean(Gini.sp)+sd(Gini.sp)*qt(0.975,m-1))
print(CI)          #The  approximate 95% conffidence interval 
cover.rate<-sum(I(Gini.sp>CI[1]&Gini.sp<CI[2]))/m
print(cover.rate)  #The result can show the CI is efficient.

## ------------------------------------------------------------------------
m<-1000
sigma1<-matrix(c(1.5,0.375,0.375,4),2,2)
p.value1=p.value2=p.value3=numeric(m)
for(i in 1:1000){
x1<-mvrnorm(20,mu=c(0,0),Sigma = sigma1)
p.value1[i]<-cor.test(x1[,1],x1[,2],method = "pearson")$p.value
p.value2[i]<-cor.test(x1[,1],x1[,2],method = "spearman")$p.value
p.value3[i]<-cor.test(x1[,1],x1[,2],method = "kendall")$p.value
}
#the power
power<-c(mean(p.value1<=0.05),mean(p.value2<=0.05),mean(p.value3<=0.05))
names(power)<-c("Pearson","Spearman","Kendall")
print(power) ##This result shows that the power of Spearman is bigger than the power of Pearson.

## ------------------------------------------------------------------------
library(boot)
library(bootstrap)
library(DAAG)

## ------------------------------------------------------------------------
set.seed(1)
attach(law)
b.cor<-function(x,i){
      y=cor(x[1,i],x[2,i])
      return(y)
}
x=t(law)
n=nrow(law)
theta.hat <- b.cor(x,1:n)
theta.jack <- numeric(n)
for(i in 1:n){
theta.jack[i] <- b.cor(x,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias=bias.jack,se=se.jack),3)

## ------------------------------------------------------------------------
attach(aircondit)
X=aircondit$hours

mle.boot<-function(X,i) mean(X[i])
de <- boot(data=X,statistic=mle.boot, R =1000)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
ci$norm[2:3]
ci$basic[4:5]
ci$percent[4:5]
ci$bca[4:5]

## ------------------------------------------------------------------------

attach(scor)
X=scor
n=dim(X)[1]
p=dim(X)[2]
eigen=numeric(p)
A<-diag(p)
theta.jac=numeric(n)
eigen=eigen(cov(X))$values
theta.hat=max(eigen)/sum(eigen)
for(i in 1:n){
   A=cov(scor[-i,])
   eigen=eigen(A)$values
   theta.jac[i]=max(eigen)/sum(eigen)
}
bias.jac <- (n-1)*(mean(theta.jac)-theta.hat)
se.jac <- sqrt((n-1)*mean((theta.jac-theta.hat)^2))
round(c(original=theta.hat,bias=bias.jack, se=se.jack),3)

## ------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- cbind(rep(0,n*(n-1)/2),rep(0,n*(n-1)/2))
k=1
for(i in 1:(n-1)){
   for(j in (i+1):n ){
         y <- magnetic[-c(i,j)]
         x <- chemical[-c(i,j)]
         J1 <- lm(y ~ x)
         yhat11 <- J1$coef[1] + J1$coef[2] * chemical[i]
         yhat12 <- J1$coef[1] + J1$coef[2] * chemical[j]
         e1[k,1]=yhat11-magnetic[i]
         e1[k,2]=yhat12-magnetic[j]
         
         J2 <- lm(y ~ x + I(x^2))
         yhat21 <- J2$coef[1] + J2$coef[2] * chemical[i] +J2$coef[3] * chemical[i]^2
         yhat22 <- J2$coef[1] + J2$coef[2] * chemical[j] +J2$coef[3] * chemical[j]^2
         e2[k,1]=yhat21-magnetic[i]
         e2[k,2]=yhat22-magnetic[j]
         J3 <- lm(log(y) ~ x)
         logyhat31 <- J3$coef[1] + J3$coef[2] * chemical[i]
         yhat31 <- exp(logyhat31)
         logyhat32 <- J3$coef[1] + J3$coef[2] * chemical[j]
         yhat32 <- exp(logyhat32)
         e3[k,1] <-  yhat31-magnetic[i]
         e3[k,2] <-  yhat32-magnetic[j] 
         J4 <- lm(log(y) ~ log(x)) 
         logyhat41 <- J4$coef[1] + J4$coef[2] * log(chemical[i])
         logyhat42 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
         yhat41 <- exp(logyhat41)
         yhat42 <- exp(logyhat42)
         e4[k,1] <-  yhat41 -magnetic[i]
         e4[k,2] <-  yhat42 -magnetic[j]
         k=k+1
    }
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## ------------------------------------------------------------------------
set.seed(1)
library(boxplotdbl)
attach(chickwts)
cvm.ts<-function(x,y){
  n=length(x);m=length(y)
  fn=ecdf(x);gn=ecdf(y)
  w2=m*n/((m+n)^2)*(sum((fn(x)-gn(x))^2)+sum((fn(y)-gn(y))^2))
  return(w2)
}
N=999
x<-sort(as.vector(weight[feed=="soybean"]))
y<-sort(as.vector(weight[feed=="linseed"]))
detach(chickwts)
z=c(x,y)
K=1:26
reps=numeric(N)
cvm0=cvm.ts(x,y)
for (i in 1:N) {
  k<-sample(K,size = 14,replace = F)
  x1<-z[k]
  y1<-z[-k]
  reps[i]<-cvm.ts(x1,y1)
}
p<-mean(c(cvm0,reps)>=cvm0)
p#0.421 does not support the alternative hypothesis that distributions differ
hist(reps,main="",freq=F,xlab="W2(p=0.421)",breaks = "scott")
points(cvm0,0,cex=1,pch=16)

## ------------------------------------------------------------------------
library(Ball)
library(energy)
library(boot)
library(RANN)

## ------------------------------------------------------------------------
m<-1e2;k<-3;p<-2;mu<-0.5;R<-999;
p.values <- matrix(0,m,3)

Tn<-function(z,ix,sizes,k){
n1<-sizes[1];n2<-sizes[2];n<-n1+n2
if(is.vector(z)) z<-data.frame(z,0)
z<-z[ix, ]
NN<-nn2(data=z,k=k+1)
block1<-NN$nn.idx[1:n1,-1] 
block2<-NN$nn.idx[(n1+1):n,-1] 
i1<-sum(block1 <= n1)
i2<-sum(block2 >n1) 
(i1+i2)/(k*n)
}

eqdist.nn<-function(z,sizes,k){
boot.obj<-boot(data=z,statistic=Tn,R=R,sim="permutation",sizes=sizes,k=k)
ts<-c(boot.obj$t0,boot.obj$t)
p.value<-mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}

## ------------------------------------------------------------------------
n1<-n2<-50;n<-n1+n2;N=c(n1,n2);mu1=mu2=0;sigma1=1;sigma2=2  
for(i in 1:m){
x<-rnorm(n1,mu1,sigma1)
y<-rnorm(n2,mu2,sigma2)
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
n1<-n2<-50;n<-n1+n2;N=c(n1,n2);mu1=1;mu2=2;sigma1=1;sigma2=2   
for(i in 1:m){
x<-rnorm(n1,mu1,sigma1)
y<-rnorm(n2,mu2,sigma2)
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
n1<-n2<-50;n<-n1+n2;N=c(n1,n2);mu1=1;mu2=2;sigma1=1;sigma2=2  
for(i in 1:m){
x<-rt(n1,1)
y<-c(rnorm(n2/2,mu1,sigma1),rnorm(n2/2,mu2,sigma2))
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
n1<-10;n2<-100;n<-n1+n2;N=c(n1,n2) 
for(i in 1:m){
x<-rnorm(n1)
y<-rnorm(n2)
z<-c(x,y)
p.values[i,1]<-eqdist.nn(z,N,k)$p.value
p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3]<-bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha<-0.1
pow<-colMeans(p.values<alpha) 
pow

## ------------------------------------------------------------------------
 f_cauchy<-function(x,u,lamada){
    return(lamada/(pi*(lamada^2+(x-u)^2)))
} #the density function of cauchy

cauchy.chain<-function(n,sigma,x0,u,lamada){
    #n:length of the chain
    #sigma:standard error of the proposal distribution
    #x0:initial value
    #u,lamada:parameter of cauchy distribution
    x<-NULL
    x[1]<-x0
    e=runif(n)
    k<-0
    for(i in 2:n){
      y<-rnorm(1,x[i-1],sigma)
      if(e[i]<=(f_cauchy(y,u,lamada)/f_cauchy(x[i-1],u,lamada)))  x[i]<-y
      else{
      x[i]<-x[i-1]
      k<-k+1
      }
    }
    return(list(x=x,k=k))
  } 
  n<-10000
  cauchy<-cauchy.chain(n,sigma=0.5,x0=0,u=0,lamada=1)
  refuse.pr<-cauchy$k/n
  refuse.pr
  #qq plot
  par(mfrow=c(1,1))
  qqplot(rcauchy(n-1000,0,1),cauchy$x[1001:n])
  qqline(cauchy$x[1001:n])
  #histplot
  hist(cauchy$x[1001:n],freq=F,main="density of cauchy",breaks=60)
  curve(f_cauchy(x,0,1),add=TRUE)


## ------------------------------------------------------------------------
  w <-0.25 #width of the uniform support set 
  n <- 5000 #length of the chain 
  burn <- 1000 #burn-in time 
  y<-c(125,18,20,34)
  x <- numeric(n) #the chain
  
  prob <- function(b, y) { 
    # computes (without the constant) the target density 
    if (b < 0 || b >= 1) return (0)
    return((1/2+b/4)^y[1] * ((1-b)/4)^y[2] * ((1-b)/4)^y[3] * (b/4)^y[4])
  }
  
  u <- runif(n) #for accept/reject step
  v <- runif(n, -w, w) #proposal distribution 
  x[1] <-0.25 
  for (i in 2:n) { 
    z <- x[i-1] + v[i] 
    if (u[i] <= prob(z,y) / prob(x[i-1],y)) 
      x[i] <-z 
    else 
      x[i] <- x[i-1] 
  }
  
   xb <- x[(burn+1):n]
   xc<-mean(xb)
   print(xc)  #estimation value of theta
   print(y/sum(y))
   print(c(1/2+xc/4,(1-xc)/4,(1-xc)/4,xc/4))


## ------------------------------------------------------------------------
set.seed(1)
Gelman.Rubin <- function(psi) {
  psi<-as.matrix(psi)        # psi[i,j] is the statistic psi(X[i,1:j])
  n<-ncol(psi)               # for chain in i-th row of X
  k<-nrow(psi)
  psi.means<-rowMeans(psi)   #row means
  B<-n*var(psi.means)        #between variance est.
  psi.w<-apply(psi,1,"var")  #within variances
  W<-mean(psi.w)             #within est.
  v.hat<-W*(n-1)/n+(B/n)     #upper variance est.
  r.hat<-v.hat/W             #G-R statistic
  return(r.hat)
}       

prob<-function(y,ob){ # computes the target density 
  if(y < 0 || y >= 1){
    return(0)
    }
  else{
    return((0.5+y/4)^ob[1]*((1-y)/4)^ob[2]*((1-y)/4)^ob[3]*(y/4)^ob[4])
    }
}

mult.chain<-function(ob,w,m,x1){
   x<-numeric(m)      #the chain
   u<-runif(m)        #for accept/reject step 
   v<-runif(m,-w,w)   #proposal distribution 
   x[1]<-x1 
   for(i in 2:m) { 
   y<-x[i-1]+v[i] 
      if(u[i]<=prob(y,ob)/prob(x[i-1],ob)){
            x[i]<-y
      }
   else{ 
            x[i]<-x[i-1]} 
   }
   return(x)
}

w<-.5               #width of the uniform support set 
ob<-c(125,18,20,34) #the observed sample
m<-10000            #length of the chain 
burn<-1000          #burn-in time 
k<-4                #number of chains to generate
x0<-c(0.5,0.9,0.1,0.75)  #choose overdispersed initial values

#generate the chains
X<-matrix(0,nrow=k,ncol=m)
for(i in 1:k){
  X[i, ]<-mult.chain(ob,w,m,x0[i])
}
    
psi<-t(apply(X,1,cumsum)) #compute diagnostic statistics
for(i in 1:nrow(psi)){
  psi[i,]<-psi[i,]/(1:ncol(psi))
}

# chain has approximately converged to the target distribution within approximately 3000 iterations 

## ------------------------------------------------------------------------
Ak<-function(a,k){     
  t1<- sqrt(a^2*(k-1)/(k-a^2))
  t2<- sqrt(a^2*k/(k+1-a^2)) 
  return(pt(t1,df=k-1)-pt(t2,df=k))
}
##Function Ak=0 means that two curves  intersect.

K=c(4:25,100,500,1000)
ipts=numeric(length(K))
for(i in 1:length(K)){
    ipts[i]=uniroot(Ak,lower=1e-4,upper=sqrt(K[i])-(1e-4),k=K[i])$root
}
ipts #output the result(Wrong)

## Check the ipts, we can find that K>=23 leads to ipts arrive at right endpoint.
## at right endpoint. Which may be wrong.
## So we should chang the right endpoint based on the curves.
k=23  ## for example
x=seq(0.01,sqrt(k)-1e-4,length=1000)
y=Ak(x,k)
plot(x,y,type="l")
## The curve shows that the null point is smaller than 3.
## For other k is the same result.
## So chang the right endpoint to 3.
wrong=which(K>22)
for (i in wrong){
  ipts[i] <- uniroot(Ak,lower =1e-4,upper=3,k=K[i])$root
}
names(ipts)=K
ipts   ## the final result

## ------------------------------------------------------------------------
library(nloptr)
f<-function(y,eta,theta){
#theta is scale parameter
#eta is the location parameter
1/(theta*3.141592653*(1+((y-eta)/theta)^2))
}

# the cauchy pdf
p_cauchy<-function(x,eta,theta,lower.tail=TRUE){
 if(lower.tail) res<-integrate(f,lower = -Inf,upper = x,rel.tol=.Machine$double.eps^0.25,theta=theta,eta=eta)
 else res<-integrate(f,lower = x,upper = Inf,rel.tol=.Machine$double.eps^0.25,theta=theta,eta=eta)
  return(res$value)
}
##compare
p_cauchy(x=0,eta = 0,theta = 1)
pcauchy(0,location = 0,scale = 1)

p_cauchy(x=1,eta =3,theta = 2,lower.tail = F )
pcauchy(1,location = 3,scale = 2,lower.tail = F)

## ----echo=FALSE----------------------------------------------------------
options(warn=-1)
dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB'),
             Frequency=c('p2','q2','r2','2pr','2qr','2pq',1),
             Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
knitr::kable(dat,format='markdown')

## ------------------------------------------------------------------------
N<-10000 #max. number of iterations 
tol<-.Machine$double.eps 
L.old<-c(.2,.35)  #initial est. for p and q 
M.value<-0
L.list<-data.frame(p=0,q=0)

mlef<-function(l,l.old,n.A=28,n.B=24,n.OO=41,n.AB=70 ){
r<-1-sum(l)
r.old<-1-sum(l.old)
n.AA<-n.A*l.old[1]^2/(l.old[1]^2+2*l.old[1]*r.old)
n.BB<-n.B*l.old[2]^2/(l.old[2]^2+2*l.old[2]*r.old)
llh<-2*n.OO*log(r)+n.AB*log(2*l[1]*l[2])+
2*n.AA*log(l[1])+2*n.BB*log(l[2])+
(n.A-n.AA)*log(2*l[1]*r)+(n.B-n.BB)*log(2*l[2]*r)
-llh
}

for(j in 1:N){
res<-optim(c(0.3,0.2),mlef,l.old=L.old)
L<-res$par
L.list[j,]<-L
M.value[j]<- -res$value
if(sum(abs(L-L.old)/L.old)<tol) break 
L.old<-L
}
L.list  #p,q
print(-M.value)  ##the max likelihood values

## ------------------------------------------------------------------------
set.seed(123)
attach(mtcars)
formulas <- list( mpg ~ disp, mpg ~ I(1 / disp), mpg ~ disp + wt, mpg ~ I(1 / disp) + wt )
#loops
output <- vector("list", length(formulas))
for (i in seq_along(formulas)){
  output[[i]]<-lm(formulas[[i]])
}
output
#lapply
lapply(formulas,lm)

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ] })

#loops
for(i in seq_along(bootstraps)){
  print(lm(mpg~disp,data =bootstraps[[i]]))
}

#lapply
lapply(bootstraps,lm,formula=mpg~disp)


## ------------------------------------------------------------------------
rsq <- function(mod) summary.lm(mod)$r.squared
#loops
for (i in seq_along(formulas)) {
 print( rsq(lm(formulas[[i]])))
  }
#lapply
lapply(lapply(formulas,lm),rsq)
#loops
for(i in seq_along(bootstraps)){
  print(rsq(lm(mpg~disp,data =bootstraps[[i]])))
}

#lapply
lapply(lapply(bootstraps,lm,formula=mpg~disp),rsq)

## ------------------------------------------------------------------------
#using anonymous function
trials <- replicate( 100, t.test(rpois(10, 10), rpois(7, 10)), simplify = FALSE )
p_value<-function(mod) mod$p.value
sapply(trials, p_value)



## ------------------------------------------------------------------------
Mvap<-function (f,n,type, ...) {  #n is the length of output of function f. type is the mode of output of function f.
  f <- match.fun(f)
  tt=Map(f, ...)
  if(type=="numeric")  y=vapply(tt,cbind,numeric(n))
  else if (type=="character") y=vapply(tt,cbind,character(n))
  else if (type=="complex") y=vapply(tt,cbind,complex(n))
  else if (type=="logical") y=vapply(tt,cbind,logical(n))
  return(y)
}
#example
kk<-function(x,w){
  y=weighted.mean(x, w, na.rm = TRUE)
  return(c(1,y))
}
xs <- replicate(5, runif(10), simplify = FALSE)
ws <- replicate(5, rpois(10, 5) + 1, simplify = FALSE)
Mvap(kk,2,"numeric",xs, ws)
is.matrix(Mvap(kk,2,"numeric",xs, ws))

## ----echo=FALSE----------------------------------------------------------
library(microbenchmark)

## ------------------------------------------------------------------------
new.chisq.test<-function(x,y){
if(!is.vector(x) && !is.vector(y))
stop("at least one of 'x' and 'y' is not a vector")
if(typeof(x)=="character" || typeof(y)=="character")
stop("at least one of 'x' and 'y' is not a numeric vector")
if(any(x<0) || anyNA(x)) 
stop("all entries of 'x' must be nonnegative and finite")
if(any(y<0) || anyNA(y)) 
stop("all entries of 'y' must be nonnegative and finite")
if((n<-sum(x))==0) 
stop("at least one entry of 'x' must be positive")
if((n<-sum(x))==0) 
stop("at least one entry of 'x' must be positive")
if(length(x)!=length(y)) 
stop("'x' and 'y' must have the same length")
DNAME<-paste(deparse(substitute(x)),"and",deparse(substitute(y)))
METHOD<-"Pearson's Chi-squared test"
x<-rbind(x,y)
nr<-as.integer(nrow(x));nc<-as.integer(ncol(x))
sr<-rowSums(x);sc<-colSums(x);n<-sum(x)
E<-outer(sr,sc,"*")/n
STATISTIC<-sum((x - E)^2/E)
names(STATISTIC)<-"X-squared"
structure(list(statistic=STATISTIC,method=METHOD,data.name=DNAME),class="htest")
}

#example
a<-c(568,327,494);b<-c(484,339,444)    
new.chisq.test(a,b)        
chisq.test(rbind(a,b))
microbenchmark(t1=new.chisq.test(a,b),t2=chisq.test(rbind(a,b)))  

## ------------------------------------------------------------------------
new.table<-function(...,dnn = list.names(...),deparse.level = 1){
    list.names <- function(...) {
        l <- as.list(substitute(list(...)))[-1L]
        nm <- names(l)
        fixup <- if (is.null(nm)) 
            seq_along(l)
        else nm == ""
        dep <- vapply(l[fixup], function(x) switch(deparse.level + 
            1, "", if (is.symbol(x)) as.character(x) else "", 
            deparse(x, nlines = 1)[1L]), "")
        if (is.null(nm)) 
            dep
        else {
            nm[fixup] <- dep
            nm
        }
    }
    args <- list(...)
    if (!length(args)) 
        stop("nothing to tabulate")
    if (length(args) == 1L && is.list(args[[1L]])) {
        args <- args[[1L]]
        if (length(dnn) != length(args)) 
            dnn <- if (!is.null(argn <- names(args))) 
                argn
            else paste(dnn[1L], seq_along(args), sep = ".")
    }
    bin <- 0L
    lens <- NULL
    dims <- integer()
    pd <- 1L
    dn <- NULL
    for (a in args) {
        if (is.null(lens)) 
            lens <- length(a)
        else if (length(a) != lens) 
            stop("all arguments must have the same length")
        fact.a <- is.factor(a)
        if (!fact.a) {
            a0 <- a
            a <- factor(a)
        }
        ll <- levels(a)
        a <- as.integer(a)
        nl <- length(ll)
        dims <- c(dims, nl)
        dn <- c(dn, list(ll))
        bin <- bin + pd * (a - 1L)
        pd <- pd * nl
    }
    names(dn) <- dnn
    bin <- bin[!is.na(bin)]
    if (length(bin)) 
        bin <- bin + 1L
    y <- array(tabulate(bin, pd), dims, dimnames = dn)
    class(y) <- "table"
    y
}

#example
a<-b<-c(1,1,2,3)            
new.table(a,b)
table(a,b)
microbenchmark(t1=new.table(a,b),t2=table(a,b))    

