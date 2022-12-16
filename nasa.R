##
## Code for analysis of nasa data from preflib 
## https://www.preflib.org/static/data/nasa/00003-00000001.toc
##
## author: daniel.henderson@newcastle.ac.uk

## read in  GPL functions
source("gpl-functions.R")

## read in data and convert to y,s format
y <- read.csv("nasa-y.txt",header=FALSE)
y <- as.matrix(y)
t <- read.csv("nasa-t.txt",header=FALSE)
t <- as.matrix(t)
s <- t.to.s(t)


## extract number of entities K, and number of comparisons n 
K <- dim(y)[2]
n <- dim(y)[1]


##################################################
## GPL model (smaller is better)

## compute summary statistics
ss <- sum.stat.gpl(K,y,s)

## beta prior hyperparameters
a <- rep(1,K)
b <- rep(1,K)

## number of iterations
its <- 10010

set.seed(100)
## run 4 chains

## initial values
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
system.time(res <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr))

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
res2 <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
res3 <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
res4 <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr)

## traceplots
ts.plot(res$theta,col=1:K)

## matrix of pointwise log-likelihoods over each datapoint
ll.res <- matrix(0,nrow=its,ncol=n)
ll2.res <- matrix(0,nrow=its,ncol=n)
ll3.res <- matrix(0,nrow=its,ncol=n)
ll4.res <- matrix(0,nrow=its,ncol=n)
for(i in 1:its){
  ll.res[i,] <- loglike.gpl(y,t,res$theta[i,],m.i.star=rep(K-1,n))
  ll2.res[i,] <- loglike.gpl(y,t,res2$theta[i,],m.i.star=rep(K-1,n))
  ll3.res[i,] <- loglike.gpl(y,t,res3$theta[i,],m.i.star=rep(K-1,n))
  ll4.res[i,] <- loglike.gpl(y,t,res4$theta[i,],m.i.star=rep(K-1,n))
}
plot.ts(ll.res)
plot.ts(ll2.res)
plot.ts(ll3.res)
plot.ts(ll4.res)

## choose burn-in

burn <- 10
its <- its-burn

## remove the first burn iterations
ll.res <- ll.res[-(1:burn),]
ll2.res <- ll2.res[-(1:burn),]
ll3.res <- ll3.res[-(1:burn),]
ll4.res <- ll4.res[-(1:burn),]
res$theta <- res$theta[-(1:burn),]
res2$theta <- res2$theta[-(1:burn),]
res3$theta <- res3$theta[-(1:burn),]
res4$theta <- res4$theta[-(1:burn),]
plot.ts(ll.res)
plot.ts(ll2.res)
plot.ts(ll3.res)
plot.ts(ll4.res)


par(mfrow=c(4,4))
for(i in 1:16){
  plot(res$theta[,i],col=1,type="l")
  lines(res2$theta[,i],col=2)
  lines(res3$theta[,i],col=3)
  lines(res4$theta[,i],col=4)
}
for(i in 1:16){
  plot(res$theta[,16+i],col=1,type="l")
  lines(res2$theta[,16+i],col=2)
  lines(res3$theta[,16+i],col=3)
  lines(res4$theta[,16+i],col=4)
}
par(mfrow=c(1,1))

## summaries
library(coda)

res.mat <- as.mcmc.list(list(as.mcmc(res$theta),as.mcmc(res2$theta),as.mcmc(res3$theta),as.mcmc(res4$theta)))

psrf <- gelman.diag(res.mat)
psrf$psrf
round(psrf$psrf,3)

## results from 1 chain
effectiveSize(res$theta)
print(min(effectiveSize(res$theta)))
which.min(effectiveSize(res$theta))
acf(res$theta[,1])
acf(res$theta[,which.min(effectiveSize(res$theta))])

pdf("nasa-mix.pdf",width=15,height=5,pointsize=16)
par(mfrow=c(1,2))
plot(res$theta[,which.min(effectiveSize(res$theta))],type="l",xlab="Iteration",ylab=expression(theta[11]),main="")
acf(res$theta[,which.min(effectiveSize(res$theta))],main="")
par(mfrow=c(1,1))
dev.off()


## results from all 4 chains
effectiveSize(res.mat)
min(effectiveSize(res.mat))
which.min(effectiveSize(res.mat))


summary(as.mcmc(res$theta))
##plot(as.mcmc(res$theta))

## low posterior correlations
##pairs(res$theta,pch=".")
crosscorr(as.mcmc(res$theta))
crosscorr.plot(as.mcmc(res$theta))


## order the entities in terms of posterior mean (largest to smallest)
theta.postmean <- apply(res$theta,2,mean)
rev(order(theta.postmean))


## check for consistency over multiple chains
cbind(rev(order(apply(res$theta,2,mean))),rev(order(apply(res2$theta,2,mean))),rev(order(apply(res3$theta,2,mean))),rev(order(apply(res4$theta,2,mean))))
## good agreement across multiple chains 

## LOOIC

library(loo)

releff <- relative_eff(exp(ll.res),chain_id=rep(1,its))

loo.res <- loo(ll.res,r_eff=releff)
loo.res$estimates
loo.res$pointwise

## WAIC 
waic(ll.res,r_eff=releff)

#############################################################
## EM algorithm

## set the seed
set.seed(2)

## maximum number of iterations
max.its <- 60

## stopping criteria for MSD - mean squared difference
stop.crit <- 10^(-16)


## initial values
theta.curr <- rbeta(K,a,b)

em.res <- em.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,max.its,stop.crit,theta.curr)

## order the entities in terms of MAP (largest to smallest)
theta.MAP <- em.res$theta[em.res$its,]
rev(order(theta.MAP))

## comparison (similar but not the same)
cbind(rev(order(theta.MAP)),rev(order(theta.postmean)))

## posterior means and 95% equi-tailed credible intervals (plus MAP values)
pdf(file="nasa-postMAP.pdf",width=8,height=5,pointsize=14)
plot(1:K,theta.postmean,ylim=c(0,0.5),pch=19,ylab=expression(theta[k]),xlab="k")
points(1:K,theta.MAP,pch=4)
for(k in 1:K){
  lines(c(k,k),quantile(res$theta[,k],c(0.025,0.975)))
}
dev.off()


#############################################################
## reverse (bigger is better) model

y.rev <- t(apply(y,1,rev))
s.rev <- s[,K]+1-t(apply(s,1,rev))
t.rev <- t(apply(t,1,rev)) ## double check

y <- y.rev
s <- s.rev
t <- t.rev

## compute summary statistics
ss <- sum.stat.gpl(K,y,s)

## beta prior hyperparameters
a <- rep(1,K)
b <- rep(1,K)

## number of iterations
its <- 10010

set.seed(200)
## run 4 chains

## initial values
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
system.time(res.rev <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr))
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
res2.rev <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr)
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
res3.rev <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr)
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
res4.rev <- gs.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,its,theta.curr,z.curr)

## traceplots
ts.plot(res.rev$theta,col=1:K)

## matrix of pointwise log-likelihoods over each datapoint
ll.res.rev <- matrix(0,nrow=its,ncol=n)
ll2.res.rev <- matrix(0,nrow=its,ncol=n)
ll3.res.rev <- matrix(0,nrow=its,ncol=n)
ll4.res.rev <- matrix(0,nrow=its,ncol=n)
for(i in 1:its){
  ll.res.rev[i,] <- loglike.gpl(y,t,res.rev$theta[i,],m.i.star=rep(K-1,n))
  ll2.res.rev[i,] <- loglike.gpl(y,t,res2.rev$theta[i,],m.i.star=rep(K-1,n))
  ll3.res.rev[i,] <- loglike.gpl(y,t,res3.rev$theta[i,],m.i.star=rep(K-1,n))
  ll4.res.rev[i,] <- loglike.gpl(y,t,res4.rev$theta[i,],m.i.star=rep(K-1,n))
}
plot.ts(ll.res.rev)
plot.ts(ll2.res.rev)
plot.ts(ll3.res.rev)
plot.ts(ll4.res.rev)


## decide on burn-in

burn <- 10
its <- its-burn

## remove first burn iterations
ll.res.rev <- ll.res.rev[-(1:burn),]
ll2.res.rev <- ll2.res.rev[-(1:burn),]
ll3.res.rev <- ll3.res.rev[-(1:burn),]
ll4.res.rev <- ll4.res.rev[-(1:burn),]
res.rev$theta <- res.rev$theta[-(1:burn),]
res2.rev$theta <- res2.rev$theta[-(1:burn),]
res3.rev$theta <- res3.rev$theta[-(1:burn),]
res4.rev$theta <- res4.rev$theta[-(1:burn),]
plot.ts(ll.res.rev)
plot.ts(ll2.res.rev)
plot.ts(ll3.res.rev)
plot.ts(ll4.res.rev)


par(mfrow=c(4,4))
for(i in 1:16){
  plot(res.rev$theta[,i],col=1,type="l")
  lines(res2.rev$theta[,i],col=2)
  lines(res3.rev$theta[,i],col=3)
  lines(res4.rev$theta[,i],col=4)
}
for(i in 1:16){
  plot(res.rev$theta[,16+i],col=1,type="l")
  lines(res2.rev$theta[,16+i],col=2)
  lines(res3.rev$theta[,16+i],col=3)
  lines(res4.rev$theta[,16+i],col=4)
}
par(mfrow=c(1,1))

## summaries
library(coda)

res.mat.rev <- as.mcmc.list(list(as.mcmc(res.rev$theta),as.mcmc(res2.rev$theta),as.mcmc(res3.rev$theta),as.mcmc(res4.rev$theta)))

psrf.rev <- gelman.diag(res.mat.rev)
psrf.rev$psrf
print(round(psrf.rev$psrf,3))

## results from 1 chain

effectiveSize(res.rev$theta)
print(min(effectiveSize(res.rev$theta)))
which.min(effectiveSize(res.rev$theta))
acf(res.rev$theta[,1])
acf(res.rev$theta[,which.min(effectiveSize(res.rev$theta))])

## results from all 4 chains
effectiveSize(res.mat.rev)
min(effectiveSize(res.mat.rev))
which.min(effectiveSize(res.mat.rev))

summary(as.mcmc(res.rev$theta))
##plot(as.mcmc(res.rev$theta))

## posterior correlations
crosscorr.plot(as.mcmc(res.rev$theta))

## order the entities in terms of posterior mean (smallest to largest)
theta.rev.postmean <- apply(res.rev$theta,2,mean)
print(round(theta.rev.postmean,3))
print(order(theta.rev.postmean))

## check for consistency over multiple chains
cbind(order(apply(res.rev$theta,2,mean)),order(apply(res2.rev$theta,2,mean)),order(apply(res3.rev$theta,2,mean)),order(apply(res4.rev$theta,2,mean)))
## very similar same

## LOOIC

library(loo)

releff.r <- relative_eff(exp(ll.res.rev),chain_id=rep(1,its))

loo.res.r <- loo(ll.res.rev,r_eff=releff.r)
print(loo.res.r$estimates)
loo.res.r$pointwise

## WAIC 
waic(ll.res.rev,r_eff=releff.r)

######################################################
## EM algorithm

## set the seed
set.seed(2)

## maximum number of iterations
max.its <- 60

## stopping criteria for MSD - mean squared difference
stop.crit <- 10^(-16)

## initial values
theta.curr <- rbeta(K,a,b)

em.rev <- em.gpl(K,y,s,v=ss$v,w=ss$w,delta=ss$delta,a,b,max.its,stop.crit,theta.curr)

theta.MAP.rev <- em.rev$theta[em.rev$its,]
print(round(theta.MAP.rev,3))

## order the entities in terms of MAP (smallest to largest)
print(order(theta.MAP.rev))

## posterior means and 95% equi-tailed credible intervals (plus true values)
pdf(file="nasa-rev-postMAP.pdf",width=8,height=5,pointsize=14)
plot(1:K,theta.rev.postmean,ylim=c(0,0.5),pch=19,ylab=expression(theta[k]),xlab="k")
points(1:K,theta.MAP.rev,pch=4)
for(k in 1:K){
  lines(c(k,k),quantile(res.rev$theta[,k],c(0.025,0.975)))
}
dev.off()


#############################################################
## Davidson-Luce

library(PlackettLuce)

## read in data and convert to y,s format
y <- read.csv("nasa-y.txt",header=FALSE)
y <- as.matrix(y)
t <- read.csv("nasa-t.txt",header=FALSE)
t <- as.matrix(t)
s <- t.to.s(t)


## extract number of entities K, and number of comparisons n 
K <- dim(y)[2]
n <- dim(y)[1]

## put in correct format (ranks)

r <- matrix(0,nrow=n,ncol=K)
for(i in 1:n){
  r[i,] <- s[i,order(y[i,])]
}

r.rankings <- as.rankings(r,input="rankings")

## count up each ranking (all unique)
w = rep(1,n)

## iterative scaling

print("Starting fitting the Davidson-Luce model at")
print(date())
dl.nasa <- PlackettLuce(r.rankings, weights = w,npseudo = 0,maxit = 1)
## time out after 30 mins





