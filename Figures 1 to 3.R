rm(list=ls())
library(mvtnorm)
library(Matrix)
graphics.off()

# This funciton calculates means, variances etc of a distribution
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

# Standard Gaussian inheritance function
H.fun <- function(z, zz, mu.a, mu.b, sigma.a) {
  mu.z <- mu.a+mu.b*z
  sigma.z2 <- sigma.a
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}

# Reproduction function -- generates offspring distribution
reproduction <- function(mum,dad,E.var){
  mum <- tapply(mum,z[,1],sum)
  dad <- tapply(dad,z[,1],sum)
  mum.moments  <- moments.fun(mum,A.vals)
  dad.moments  <- moments.fun(dad,A.vals)
  mid.point    <- (mum.moments[1]+dad.moments[1])/2
  mum <- rep(mum,times=n)
  dad <- rep(dad,each=n)
  midpts <- mum*dad
  mum.index <- rep(1:n,times=n)
  dad.index <- rep(1:n,each=n)
  midpt.index <- (mum.index+dad.index)/2
  added <- tapply(midpts,midpt.index,sum)
  indexer <- unique(midpt.index)
  odds <- seq.int(1,2*n-1,length.out=n)
  evens <- seq.int(2,2*n-2,length.out=n-1)
  midpts.scaled <- added[odds] + c(0,added[evens]/2) + c(added[evens]/2,0)
  tot.va       <- (mum.moments[2]+(mid.point-mum.moments[1])^2+dad.moments[2]+(mid.point-dad.moments[1])^2)/2
  mu.intercept <- 0
  mu.slope     <- 1
  va.intercept <- 0.5*tot.va
  offs.mat     <- t(outer(A.vals,A.vals,H.fun,mu.intercept,mu.slope,va.intercept))
  offs.mat     <- offs.mat/matrix(as.vector(apply(offs.mat,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  offs         <- offs.mat%*%midpts.scaled
  p1           <- rep(offs,each=n)/n
  p2           <- dnorm(E.vals,0,sqrt(E.var))
  p2           <- p2/sum(p2)*n
  p2           <- rep(p2,times=n)
  inhe         <- (p1 * p2)
  return(inhe)
}

runner.if <- function(E.var,A.var,thresh){
  Z <- array(NA,c(n*n,n.gens))
  sigma <- matrix(c(A.var,0,0,E.var),2,2)
  Z[,1] <- dmvnorm(z,mean=mu,sigma=sigma)
  Z[,1] <- Z[,1]/sum(Z[,1])
  for (i in 2:n.gens){
    females <- males <- Z[,i-1]
    mean.males <- moments.fun(males,Z.vals)[1]
    males <- ifelse(Z.vals>mean.males,thresh*Z[,i-1],Z[,i-1])
    males <- males/sum(males)
    Z[,i] <- reproduction(females,males,E.var)
  }
  return(Z)
}


# set up global parameters and choose initial values
# number of classes in the matrix
n <- 500
# values of breeding value distribution
A.vals <- seq(30,100,length.out=n)
# create the values of a in the bivariate distribution
a1 <- rep(A.vals,each=n)
# now repeat but for E
E.vals <- seq(-10,10,length.out=n)
e1 <- rep(E.vals,times=n)
# create the bivariate distribution of a and e
z <- cbind(a1,e1)
Z.vals <- z[,1]+z[,2]
# number of generations for the simulation
n.gens <- 100
# define arrays into which to place results
# start with the initial distributions
mu <- c(70,0)

# Figure 1
res.1 <- runner.if(2,3,0.75)
res.2 <- runner.if(2,3,0.5)
res.3 <- runner.if(2,3,0.25)
res.4 <- runner.if(2,3,0)

res.5 <- runner.if(0.01,4.99,0.25) # 0.99
res.6 <- runner.if(1.25,3.75,0.25) # 0.75
res.7 <- runner.if(2.5,2.5,0.25) # 0.5
res.8 <- runner.if(3.75,1.25,0.25) # 0.25

pdf('~/dropbox/JWM Bighorn sheep/figure1.pdf',useDingbats=FALSE)
#quartz()
par(mfrow=c(2,2))
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],Z.vals)[1]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic mean',main='(A) Dynamics of the\nmean phenotype',ylim=c(50,70),bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],Z.vals)[1]
lines(1:n.gens,res)
abline(h=70-1.96*sqrt(13),col='gray50',lty=3)
# Figure 1
text(96,19.5+44.5,'0.25')
text(52,19.5+44.5,'0.5')
text(96,17.6+43.5,'0.75')
text(19,19.5+44.5,'1')

abline(h=50,col='gray50',lty=2)
abline(70,-4,col='gray50',lty=2)
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],Z.vals)[2]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic variance',main='(B) Dynamics of the\nphenotypic variance',bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],Z.vals)[2]
lines(1:n.gens,res)
# Figure 1
text(95,4.8,'0.25')
text(95,3.3,'0.5')
text(95,2.6,'0.75')
text(20,2.37,'1')
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],z[,1])[2]/moments.fun(res.4[,i],Z.vals)[2]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Heritability',main='(C) Dynamics of\nthe heritability',bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],z[,1])[2]/moments.fun(res.3[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],z[,1])[2]/moments.fun(res.2[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],z[,1])[2]/moments.fun(res.1[,i],Z.vals)[2]
lines(1:n.gens,res)
# Figure 1
# NEED CHANGING......

text(95,0.575,'0.25')
text(95,0.375,'0.5')
text(95,0.2,'0.75')
text(95,0.1,'1')





res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.5[,i],Z.vals)[1]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic mean',main='(D) Dynamics of the\nmean phenotype',bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.6[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.7[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.8[,i],Z.vals)[1]
lines(1:n.gens,res)
text(95,61,'0.25')
text(20,61,'0.99')
abline(h=70-1.96*sqrt(13),col='gray50',lty=3)
dev.off()













#var.A1 <- var.Z1 <- means1 <- rep(NA,n.gens)
#for (i in 1:n.gens) means1[i] <- moments.fun(res.1[,i],z[,1])[1]
#for (i in 1:n.gens) var.A1[i] <- moments.fun(res.1[,i],z[,1])[2]
#for (i in 1:n.gens) var.Z1[i] <- moments.fun(res.1[,i],Z.vals)[2]
#h21 <- var.A1/var.Z1
#delta.z1 <- c(means1[-1]-means1[-n.gens],NA)
#s1 <- delta.z1/h21

#var.A2 <- var.Z2 <- means2 <- rep(NA,n.gens)
#for (i in 1:n.gens) means2[i] <- moments.fun(res.2[,i],z[,1])[1]
#for (i in 1:n.gens) var.A2[i] <- moments.fun(res.2[,i],z[,1])[2]
#for (i in 1:n.gens) var.Z2[i] <- moments.fun(res.2[,i],Z.vals)[2]
#h22 <- var.A2/var.Z2
#delta.z2 <- c(means2[-1]-means2[-n.gens],NA)
#s2 <- delta.z2/h22

#var.A3 <- var.Z3 <- means3 <- rep(NA,n.gens)
#for (i in 1:n.gens) means3[i] <- moments.fun(res.3[,i],z[,1])[1]
#for (i in 1:n.gens) var.A3[i] <- moments.fun(res.3[,i],z[,1])[2]
#for (i in 1:n.gens) var.Z3[i] <- moments.fun(res.3[,i],Z.vals)[2]
#h23 <- var.A3/var.Z3
#delta.z3 <- c(means3[-1]-means3[-n.gens],NA)
#s3 <- delta.z3/h23

#var.A4 <- var.Z4 <- means4 <- rep(NA,n.gens)
#for (i in 1:n.gens) means4[i] <- moments.fun(res.4[,i],z[,1])[1]
#for (i in 1:n.gens) var.A4[i] <- moments.fun(res.4[,i],z[,1])[2]
#for (i in 1:n.gens) var.Z4[i] <- moments.fun(res.4[,i],Z.vals)[2]
#h24 <- var.A4/var.Z4
#delta.z4 <- c(means4[-1]-means4[-n.gens],NA)
#s4 <- delta.z4/h24


#quartz()
#par(mfrow=c(2,2))
#plot(1:100,s1,ylim=c(-1,0),xlab='Generation',ylab='Selection differential',type='l')
#lines(1:100,s2,col='black')
#lines(1:100,s3,col='black')
#lines(1:100,s4,col='black')
#text(95,-0.1,'0.25')
#text(95,-0.21,'0.5')
#text(95,-0.34,'0.75')
#text(95,-0.55,'1')

#plot(1:100,var.A1,col='black',pch=16,ylim=c(0,5),xlab='Generation',ylab='Additive genetic variance',type='l')
#lines(1:100,var.A2,col='black')
#lines(1:100,var.A3,col='black')
#lines(1:100,var.A4,col='black')

#text(95,4.2,'0.25')
#text(95,1.8,'0.5')
#text(95,0.75,'0.75')
#text(13,0.75,'1')

#Figure 2

res.1 <- runner.if(0.1,3,0.75)
res.2 <- runner.if(0.1,3,0.5)
res.3 <- runner.if(0.1,3,0.25)
res.4 <- runner.if(0.1,3,0)

pdf('~/dropbox/JWM Bighorn sheep/figure2.pdf',useDingbats=FALSE)
#quartz()
par(mfrow=c(2,2))
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],Z.vals)[1]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic mean',main='(A) Dynamics of the\nmean phenotype',ylim=c(50,70),bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],Z.vals)[1]
lines(1:n.gens,res)
abline(h=70-1.96*sqrt(13),col='gray50',lty=3)
text(50,67,'0.25')
text(96,56,'0.5')
text(40,61,'0.75')
text(10,65,'1')


abline(h=50,col='gray50',lty=2)
abline(70,-4,col='gray50',lty=2)
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],Z.vals)[2]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic variance',main='(B) Dynamics of the\nphenotypic variance',bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],Z.vals)[2]
lines(1:n.gens,res)
text(95,2.5,'0.25')
text(95,0.7,'0.5')
text(30,0.7,'0.75')
text(8,0.4,'1')
dev.off()


# Figure 3
res.1 <- runner.if(2,5,0.75)
res.2 <- runner.if(2,5,0.5)
res.3 <- runner.if(2,5,0.25)
res.4 <- runner.if(2,5,0)

pdf('~/dropbox/JWM Bighorn sheep/figure3.pdf',useDingbats=FALSE)
#quartz()
par(mfrow=c(2,2))
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],Z.vals)[1]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic mean',main='(A) Dynamics of the\nmean phenotype',ylim=c(50,70),bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],Z.vals)[1]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],Z.vals)[1]
lines(1:n.gens,res)
abline(h=70-1.96*sqrt(13),col='gray50',lty=3)
text(96,62,'0.25')
text(96,58.5,'1')
text(96,56.5,'0.75')
text(96,54,'0.5')

abline(h=50,col='gray50',lty=2)
abline(70,-4,col='gray50',lty=2)
res <- rep(NA,n.gens)
for (i in 1:n.gens) res[i] <- moments.fun(res.4[,i],Z.vals)[2]
plot(1:n.gens,res,type='l',xlab='Generations',ylab='Phenotypic variance',main='(B) Dynamics of the\nphenotypic variance',bty='L')
for (i in 1:n.gens) res[i] <- moments.fun(res.3[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.2[,i],Z.vals)[2]
lines(1:n.gens,res)
for (i in 1:n.gens) res[i] <- moments.fun(res.1[,i],Z.vals)[2]
lines(1:n.gens,res)
text(95,6,'0.25')
text(95,3.6,'0.5')
text(95,2.7,'0.75')
text(15,2.5,'1')

dev.off()
