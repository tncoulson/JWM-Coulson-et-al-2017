rm(list=ls())
graphics.off()
library(mvtnorm)
library(Matrix)
# inheritance function H(z'|z,t)  zz = z' for R code

# number of classes in the matrix
n <- 100
# values of breeding value distribution
x <- seq(-5,17,length.out=n)
z <- x
# create the combinations of x and y
x1 <- rep(x,each=n)
x2 <- rep(x,times=n)

X.vals <- cbind(x1,x2)

moments.fun <- function(x,z){
  mu <- sum(z*x)/sum(x)
  va <- sum(z*z*x)/sum(x)-mu^2
  sk <- (sum(z*z*z*x)/sum(x)-3*mu*va-mu^3)/(sqrt(va)^3)
  m4 <- sum(z*z*z*z*x)/sum(x)
  m3 <- sum(z*z*z*x)/sum(x)
  m2 <- sum(z*z*x)/sum(x)
  ku <- (m4-4*mu*m3+6*(mu^2)*m2-3*(mu^4))/(va^2)-3
  res <- c(mu,va,sk,ku,m2,m3,m4)
  names(res) <- c('Mean','Variance','Skew','Kurtosis','Moment 2','Moment 3','Moment 4')
  return(res)
}

covar <- function(x,y,n) sum(x*y*n)/sum(n)-((sum(x*n)/sum(n))*(sum(y*n)/sum(n)))

plotter <- function(mu.a1,mu.e,sigma.a1,sigma.e){
  mu.z1 <- mu.a1+mu.e
  sigma.z1 <- sigma.a1+sigma.e

  # set up global parameters and choose initial values
  # define arrays into which to place results
  Z <- A <- array(NA,c(n*n))
  # start with the initial distributions

  # distribution of breeding values, environmental component and phenotype at time 1
  A <- dmvnorm(X.vals,mean=mu.a1,sigma=sigma.a1)
  # scale to sum to unity
  A <- A/sum(A)
  Z <- dmvnorm(X.vals,mean=mu.z1,sigma=sigma.z1)
  # scale to sum to unity
  Z <- Z/sum(Z)

  w <- 0.3 + 0.1*X.vals[,1]+0.1*X.vals[,2]
  z.post.sel <- Z*w
  z.post.sel <- z.post.sel/sum(z.post.sel)
  a.post.sel <- A*w
  a.post.sel <- a.post.sel/sum(a.post.sel)

  a1.moms <- moments.fun(a.post.sel,X.vals[,1])
  a2.moms <- moments.fun(a.post.sel,X.vals[,2])
  as.covs <- covar(X.vals[,1],X.vals[,2],a.post.sel)

  mu.as <- c(a1.moms[1],a2.moms[1])
  sigma.as <- c(a1.moms[2],rep(as.covs,2),a2.moms[2])
  sigma.as <- matrix(sigma.as,2,2)

  mu.z2 <- mu.as+mu.e
  sigma.z2 <- sigma.as + sigma.e
  Z.off <- dmvnorm(X.vals,mean=mu.z2,sigma=sigma.z2)
  m1 <- moments.fun(Z,X.vals[,1])
  m2 <- moments.fun(Z,X.vals[,2])
  m1o <- moments.fun(Z.off,X.vals[,1])
  m2o <- moments.fun(Z.off,X.vals[,2])

  m1s <- moments.fun(z.post.sel,X.vals[,1])
  m2s <- moments.fun(z.post.sel,X.vals[,2])

  um1s <- (m1s[1]-m1[1])*sigma.a1[1,1]/sigma.z1[1,1]
  um2s <- (m2s[1]-m2[1])*sigma.a1[2,2]/sigma.z1[2,2]

  x <- c(m1[1],m1o[1],m1s[1],um1s)
  y <- c(m2[1],m2o[1],m2s[1],um2s)
  return(list(x,y))
}

# no environmental variance

mu.a1 <- c(6,6)
mu.e <- c(0,0)
sigma.a1 <- matrix(c(2,0,0,1),2,2)
sigma.e <- matrix(c(0.01,0,0,0.01),2,2)
evo1 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,-1.41,-1.41,1),2,2)
evo2 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,1.41,1.41,1),2,2)
evo3 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)

pdf('~/dropbox/JWM Bighorn sheep/figure4.pdf',useDingbats=FALSE)
par(mfrow=c(2,2))
evo1 <- matrix(unlist(evo1),2,4,byrow=TRUE)
plot(evo1[1,1:2],evo1[2,1:2],type='l',xlim=c(5.95,6.4),ylim=c(5.95,6.4),lty=2,xlab='Phenotype 1',ylab='Phenotype 2',bty='L')
evo2 <- matrix(unlist(evo2),2,4,byrow=TRUE)
lines(evo2[1,1:2],evo2[2,1:2],col='red',lty=2)
evo3 <- matrix(unlist(evo3),2,4,byrow=TRUE)
lines(evo3[1,1:2],evo3[2,1:2],col='blue',lty=2)

lines(evo1[1,c(1,3)],evo1[2,c(1,3)])
lines(evo2[1,c(1,3)],evo2[2,c(1,3)],col='red')
lines(evo3[1,c(1,3)],evo3[2,c(1,3)],col='blue')

abline(v=evo1[1,1]+evo1[1,4],lty=4)
abline(v=evo2[1,1]+evo2[1,4],lty=4,col='red')
abline(v=evo3[1,1]+evo3[1,4],lty=4,col='blue')

abline(h=evo1[2,1]+evo1[2,4],lty=4)
abline(h=evo2[2,1]+evo2[2,4],lty=4,col='red')
abline(h=evo3[2,1]+evo3[2,4],lty=4,col='blue')
text(5.98,6.37,'(A)')
#leg.txt <- c('Cov(A1,A2)=0','Cov(A1,A2)>0','Cov(A1,A2)<0')
#leg.col <- c('black','blue','red')
#leg.lty <- c(1,1,1)
#legend('topright',legend=leg.txt,col=leg.col,lty=leg.lty,bty='n')

# intermediate environmental variance
mu.a1 <- c(6,6)
mu.e <- c(0,0)
sigma.a1 <- matrix(c(2,0,0,1),2,2)
sigma.e <- matrix(c(2,0,0,2),2,2)
evo1 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,-1.41,-1.41,1),2,2)
evo2 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,1.41,1.41,1),2,2)
evo3 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)

evo1 <- matrix(unlist(evo1),2,4,byrow=TRUE)
plot(evo1[1,1:2],evo1[2,1:2],type='l',xlim=c(5.95,6.4),ylim=c(5.95,6.4),lty=2,xlab='Phenotype 1',ylab='Phenotype 2',bty='L')
evo2 <- matrix(unlist(evo2),2,4,byrow=TRUE)
lines(evo2[1,1:2],evo2[2,1:2],col='red',lty=2)
evo3 <- matrix(unlist(evo3),2,4,byrow=TRUE)
lines(evo3[1,1:2],evo3[2,1:2],col='blue',lty=2)

lines(evo1[1,c(1,3)],evo1[2,c(1,3)])
lines(evo2[1,c(1,3)],evo2[2,c(1,3)],col='red')
lines(evo3[1,c(1,3)],evo3[2,c(1,3)],col='blue')

abline(v=evo1[1,1]+evo1[1,4],lty=4)
abline(v=evo2[1,1]+evo2[1,4],lty=4,col='red')
abline(v=evo3[1,1]+evo3[1,4],lty=4,col='blue')

abline(h=evo1[2,1]+evo1[2,4],lty=4)
abline(h=evo2[2,1]+evo2[2,4],lty=4,col='red')
abline(h=evo3[2,1]+evo3[2,4],lty=4,col='blue')
text(5.98,6.37,'(B)')




# intermediate environmental variance
mu.a1 <- c(6,6)
mu.e <- c(0,0)
sigma.a1 <- matrix(c(2,0,0,1),2,2)
sigma.e <- matrix(c(2,1,1,2),2,2)
evo1 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,-1.41,-1.41,1),2,2)
evo2 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,1.41,1.41,1),2,2)
evo3 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)

evo1 <- matrix(unlist(evo1),2,4,byrow=TRUE)
plot(evo1[1,1:2],evo1[2,1:2],type='l',xlim=c(5.95,6.4),ylim=c(5.95,6.4),lty=2,xlab='Phenotype 1',ylab='Phenotype 2',bty='L')
evo2 <- matrix(unlist(evo2),2,4,byrow=TRUE)
lines(evo2[1,1:2],evo2[2,1:2],col='red',lty=2)
evo3 <- matrix(unlist(evo3),2,4,byrow=TRUE)
lines(evo3[1,1:2],evo3[2,1:2],col='blue',lty=2)

lines(evo1[1,c(1,3)],evo1[2,c(1,3)])
lines(evo2[1,c(1,3)],evo2[2,c(1,3)],col='red')
lines(evo3[1,c(1,3)],evo3[2,c(1,3)],col='blue')

abline(v=evo1[1,1]+evo1[1,4],lty=4)
abline(v=evo2[1,1]+evo2[1,4],lty=4,col='red')
abline(v=evo3[1,1]+evo3[1,4],lty=4,col='blue')

abline(h=evo1[2,1]+evo1[2,4],lty=4)
abline(h=evo2[2,1]+evo2[2,4],lty=4,col='red')
abline(h=evo3[2,1]+evo3[2,4],lty=4,col='blue')
text(5.98,6.37,'(C)')




# intermediate environmental variance
mu.a1 <- c(6,6)
mu.e <- c(0,0)
sigma.a1 <- matrix(c(2,0,0,1),2,2)
sigma.e <- matrix(c(2,-1,-1,2),2,2)
evo1 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,-1.41,-1.41,1),2,2)
evo2 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)
sigma.a1 <- matrix(c(2,1.41,1.41,1),2,2)
evo3 <- plotter(mu.a1,mu.e,sigma.a1,sigma.e)

evo1 <- matrix(unlist(evo1),2,4,byrow=TRUE)
plot(evo1[1,1:2],evo1[2,1:2],type='l',xlim=c(5.95,6.4),ylim=c(5.95,6.4),lty=2,xlab='Phenotype 1',ylab='Phenotype 2',bty='L')
evo2 <- matrix(unlist(evo2),2,4,byrow=TRUE)
lines(evo2[1,1:2],evo2[2,1:2],col='red',lty=2)
evo3 <- matrix(unlist(evo3),2,4,byrow=TRUE)
lines(evo3[1,1:2],evo3[2,1:2],col='blue',lty=2)

lines(evo1[1,c(1,3)],evo1[2,c(1,3)])
lines(evo2[1,c(1,3)],evo2[2,c(1,3)],col='red')
lines(evo3[1,c(1,3)],evo3[2,c(1,3)],col='blue')

abline(v=evo1[1,1]+evo1[1,4],lty=4)
abline(v=evo2[1,1]+evo2[1,4],lty=4,col='red')
abline(v=evo3[1,1]+evo3[1,4],lty=4,col='blue')

abline(h=evo1[2,1]+evo1[2,4],lty=4)
abline(h=evo2[2,1]+evo2[2,4],lty=4,col='red')
abline(h=evo3[2,1]+evo3[2,4],lty=4,col='blue')
text(5.98,6.37,'(D)')
dev.off()