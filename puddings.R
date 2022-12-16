## Analysis of puddings data from Davidson (1970)

## author: daniel.henderson@newcastle.ac.uk

## read in useful functions
source("gpl-functions.R")

## read in the data 
puddings <- read.table("puddings.txt")
## columns are i j n_ij  w_ij  w_ji  t_ij

n.row <- dim(puddings)[1]

## convert data to y and s
n <- sum(puddings[,3])
K <- 6

y <- matrix(0,nrow=n,ncol=2)
s <- matrix(0,nrow=n,ncol=2)

k <- 0
for(m in 1:n.row){
  i <- puddings[m,1]
  j <- puddings[m,2]
  for(ell in 1:puddings[m,4]){
    k <- k+1
    y[k,1] <- i
    y[k,2] <- j
    s[k,1] <- 1
    s[k,2] <- 2
  }
  for(ell in 1:puddings[m,5]){
    k <- k+1
    y[k,1] <- j
    y[k,2] <- i
    s[k,1] <- 1
    s[k,2] <- 2
  }
  for(ell in 1:puddings[m,6]){
    k <- k+1
    y[k,1] <- i
    y[k,2] <- j
    s[k,1] <- 1
    s[k,2] <- 1
  }
}      

t <- s.to.t(s)

##################################################
## GPL model (smaller is better)

## compute summary statistics
ss <- sum.stat.gpl(K,y,s)

## beta prior hyperparameters
a <- rep(1,K)
b <- rep(1,K)


## number of iterations
its <- 10010

theta.curr <- rbeta(K,1,1)
z.curr <- matrix(0,nrow=n,ncol=K) 

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
  ll.res[i,] <- loglike.gpl(y,t,res$theta[i,],m.i.star=rep(2-1,n))
  ll2.res[i,] <- loglike.gpl(y,t,res2$theta[i,],m.i.star=rep(2-1,n))
  ll3.res[i,] <- loglike.gpl(y,t,res3$theta[i,],m.i.star=rep(2-1,n))
  ll4.res[i,] <- loglike.gpl(y,t,res4$theta[i,],m.i.star=rep(2-1,n))
}
par(mfrow=c(2,2))
plot.ts(apply(ll.res,1,sum))
plot.ts(apply(ll2.res,1,sum))
plot.ts(apply(ll3.res,1,sum))
plot.ts(apply(ll4.res,1,sum))
par(mfrow=c(1,1))

## decide on amount of burn-in

burn <- 10
its <- its-burn


## chop off the first burn iterations
ll.res <- ll.res[-(1:burn),]
ll2.res <- ll2.res[-(1:burn),]
ll3.res <- ll3.res[-(1:burn),]
ll4.res <- ll4.res[-(1:burn),]
res$theta <- res$theta[-(1:burn),]
res2$theta <- res2$theta[-(1:burn),]
res3$theta <- res3$theta[-(1:burn),]
res4$theta <- res4$theta[-(1:burn),]
par(mfrow=c(2,2))
plot.ts(apply(ll.res,1,sum))
plot.ts(apply(ll2.res,1,sum))
plot.ts(apply(ll3.res,1,sum))
plot.ts(apply(ll4.res,1,sum))
par(mfrow=c(1,1))


par(mfrow=c(3,2))
for(i in 1:K){
  plot(res$theta[,i],col=1,type="l")
  lines(res2$theta[,i],col=2)
  lines(res3$theta[,i],col=3)
  lines(res4$theta[,i],col=4)
}
par(mfrow=c(1,1))


## summaries
library(coda)

res.mat <- as.mcmc.list(list(as.mcmc(res$theta),as.mcmc(res2$theta),as.mcmc(res3$theta),as.mcmc(res4$theta)))

psrf <- gelman.diag(res.mat)
psrf$psrf
print(round(psrf$psrf,3))


## results from 1 chain
effectiveSize(res$theta)
print(min(effectiveSize(res$theta)))
which.min(effectiveSize(res$theta))
acf(res$theta[,1])
acf(res$theta[,which.min(effectiveSize(res$theta))])

## results from all 4 chains
effectiveSize(res.mat)
min(effectiveSize(res.mat))
which.min(effectiveSize(res.mat))


summary(as.mcmc(res$theta))
##plot(as.mcmc(res$theta))

pairs(res$theta,pch=".")
crosscorr(as.mcmc(res$theta))
crosscorr.plot(as.mcmc(res$theta))


plot(1:K,apply(res$theta,2,mean),ylim=c(0.3,0.6),pch=19,ylab=expression(theta[k]),xlab="k")
for(k in 1:K){
  lines(c(k,k),quantile(res$theta[,k],c(0.025,0.975)))
}


##################################################################
## Compare to the alternative negative-binomial based Gibbs sampler

## compute n.mat - number of comparisons involving entities i and j
n.mat <- matrix(0,nrow=K,ncol=K)
for(i in 1:n){
  n.mat[y[i,1],y[i,2]] <- n.mat[y[i,1],y[i,2]]+1
}
n.mat <- n.mat+t(n.mat)

## Gibbs sampler
a <- rep(1,K)
b <- rep(1,K)

## number of iterations
its <- 10010

set.seed(100)
## run 4 chains

theta.curr <- rbeta(K,a,b) 
z.curr <- matrix(0,nrow=K,ncol=K)
system.time(res.nb <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr))

theta.curr <- rbeta(K,a,b) 
z.curr <- matrix(0,nrow=K,ncol=K)
res2.nb <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b) 
z.curr <- matrix(0,nrow=K,ncol=K)
res3.nb <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b) 
z.curr <- matrix(0,nrow=K,ncol=K)
res4.nb <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr)

## about 10 times quicker

## brief check posteriors are the same
par(mfrow=c(2,3))
for(i in 1:K){
  plot(density(res$theta[,i],adjust=2))
  lines(density(res.nb$theta[,i],adjust=2),col=2)
}
par(mfrow=c(1,1))


## matrix of pointwise log-likelihoods over each datapoint
ll.nb.res <- matrix(0,nrow=its,ncol=n)
ll2.nb.res <- matrix(0,nrow=its,ncol=n)
ll3.nb.res <- matrix(0,nrow=its,ncol=n)
ll4.nb.res <- matrix(0,nrow=its,ncol=n)
for(i in 1:its){
  ll.nb.res[i,] <- loglike.gpl(y,t,res.nb$theta[i,],m.i.star=rep(2-1,n))
  ll2.nb.res[i,] <- loglike.gpl(y,t,res2.nb$theta[i,],m.i.star=rep(2-1,n))
  ll3.nb.res[i,] <- loglike.gpl(y,t,res3.nb$theta[i,],m.i.star=rep(2-1,n))
  ll4.nb.res[i,] <- loglike.gpl(y,t,res4.nb$theta[i,],m.i.star=rep(2-1,n))
}
par(mfrow=c(2,2))
plot.ts(apply(ll.nb.res,1,sum))
plot.ts(apply(ll2.nb.res,1,sum))
plot.ts(apply(ll3.nb.res,1,sum))
plot.ts(apply(ll4.nb.res,1,sum))
par(mfrow=c(1,1))

## decide on burn-in
burn <- 10
its <- its-burn


## chop off the first burn iterations
ll.nb.res <- ll.nb.res[-(1:burn),]
ll2.nb.res <- ll2.nb.res[-(1:burn),]
ll3.nb.res <- ll3.nb.res[-(1:burn),]
ll4.nb.res <- ll4.nb.res[-(1:burn),]
res.nb$theta <- res.nb$theta[-(1:burn),]
res2.nb$theta <- res2.nb$theta[-(1:burn),]
res3.nb$theta <- res3.nb$theta[-(1:burn),]
res4.nb$theta <- res4.nb$theta[-(1:burn),]
par(mfrow=c(2,2))
plot.ts(apply(ll.nb.res,1,sum))
plot.ts(apply(ll2.nb.res,1,sum))
plot.ts(apply(ll3.nb.res,1,sum))
plot.ts(apply(ll4.nb.res,1,sum))
par(mfrow=c(1,1))


par(mfrow=c(3,2))
for(i in 1:K){
  plot(res.nb$theta[,i],col=1,type="l")
  lines(res2.nb$theta[,i],col=2)
  lines(res3.nb$theta[,i],col=3)
  lines(res4.nb$theta[,i],col=4)
}
par(mfrow=c(1,1))

## summaries
library(coda)

res.nb.mat <- as.mcmc.list(list(as.mcmc(res.nb$theta),as.mcmc(res2.nb$theta),as.mcmc(res3.nb$theta),as.mcmc(res4.nb$theta)))

psrf <- gelman.diag(res.nb.mat)
psrf$psrf
print(round(psrf$psrf,3))

## results from 1 chain
effectiveSize(res.nb$theta)
print(min(effectiveSize(res.nb$theta)))
which.min(effectiveSize(res.nb$theta))
acf(res.nb$theta[,1])
acf(res.nb$theta[,which.min(effectiveSize(res.nb$theta))])

pdf("pud-mix.pdf",width=15,height=5,pointsize=16)
par(mfrow=c(1,2))
plot(res.nb$theta[,which.min(effectiveSize(res.nb$theta))],type="l",xlab="Iteration",ylab=expression(theta[6]),main="")
acf(res.nb$theta[,which.min(effectiveSize(res.nb$theta))],main="")
par(mfrow=c(1,1))
dev.off()

## mixing/acf is similar to that from standard Gibbs sampler

## results from all 4 chains
effectiveSize(res.nb.mat)
min(effectiveSize(res.nb.mat))
which.min(effectiveSize(res.nb.mat))


## focus on single chain

summary(as.mcmc(res.nb$theta))
##plot(as.mcmc(res.nb$theta))

pairs(res.nb$theta,pch=".")

crosscorr(as.mcmc(res.nb$theta))
crosscorr.plot(as.mcmc(res.nb$theta))


## posterior means
theta.postmean.nb <- apply(res.nb$theta,2,mean)
print(round(theta.postmean.nb,3))

## order the entities in terms of posterior mean (largest to smallest)
print(rev(order(theta.postmean.nb)))


## LOOIC

library(loo)

releff.nb <- relative_eff(exp(ll.nb.res),chain_id=rep(1,its))

loo.nb.res <- loo(ll.nb.res,r_eff=releff.nb)
print(loo.nb.res$estimates)
loo.nb.res$pointwise

## WAIC
waic(ll.nb.res,r_eff=releff.nb)


###############################################
## EM algorithm

## set the seed
set.seed(2)

## maximum number of iterations
max.its <- 60

## stopping criteria for MSD - mean squared difference
stop.crit <- 10^(-16)

## initial values
theta.curr <- rbeta(K,a,b)

em.nb.res <- em.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,max.its,stop.crit,theta.curr)

theta.MAP.nb <- em.nb.res$theta[em.nb.res$its,]
print(round(theta.MAP.nb,3))

## order the entities in terms of MAP (largest to smallest)
print(rev(order(theta.MAP.nb)))


## posterior means and 95% equi-tailed credible intervals (plus MAP values)
pdf(file="pud-postMAP.pdf",width=8,height=5,pointsize=14)
plot(1:K,theta.postmean.nb,ylim=c(0.3,0.6),pch=19,ylab=expression(theta[k]),xlab="k")
points(1:K,theta.MAP.nb,pch=4)
for(k in 1:K){
  lines(c(k,k),quantile(res.nb$theta[,k],c(0.025,0.975)))
}
dev.off()

#######################################################################
##
## Reverse model (bigger is better) model
## 


## reverse the data 
y.rev <- t(apply(y,1,rev))
s.rev <- s[,2]+1-t(apply(s,1,rev))
t.rev <- matrix(apply(t,1,rev)) 

y <- y.rev
s <- s.rev
t <- t.rev

## compute summary statistics
ss <- sum.stat.gpl(K,y,s)

## compute n.mat - number of comparisons involving entities i and j
n.mat <- matrix(0,nrow=K,ncol=K)
for(i in 1:n){
  n.mat[y[i,1],y[i,2]] <- n.mat[y[i,1],y[i,2]]+1
}
n.mat <- n.mat+t(n.mat)


## beta prior hyperparameters
a <- rep(1,K)
b <- rep(1,K)

## number of iterations
its <- 10010

set.seed(200)
## run 4 chains
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=K,ncol=K) 
system.time(res.rev <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr))

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=K,ncol=K) 
res2.rev <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=K,ncol=K) 
res3.rev <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=K,ncol=K) 
res4.rev <- gs.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,its,theta.curr,z.curr)

## traceplots
ts.plot(res.rev$theta,col=1:K)

## matrix of pointwise log-likelihoods over each datapoint
ll.res.rev <- matrix(0,nrow=its,ncol=n)
ll2.res.rev <- matrix(0,nrow=its,ncol=n)
ll3.res.rev <- matrix(0,nrow=its,ncol=n)
ll4.res.rev <- matrix(0,nrow=its,ncol=n)
for(i in 1:its){
  ll.res.rev[i,] <- loglike.gpl(y,t,res.rev$theta[i,],m.i.star=rep(2-1,n))
  ll2.res.rev[i,] <- loglike.gpl(y,t,res2.rev$theta[i,],m.i.star=rep(2-1,n))
  ll3.res.rev[i,] <- loglike.gpl(y,t,res3.rev$theta[i,],m.i.star=rep(2-1,n))
  ll4.res.rev[i,] <- loglike.gpl(y,t,res4.rev$theta[i,],m.i.star=rep(2-1,n))
}
par(mfrow=c(2,2))
plot.ts(apply(ll.res.rev,1,sum))
plot.ts(apply(ll2.res.rev,1,sum))
plot.ts(apply(ll3.res.rev,1,sum))
plot.ts(apply(ll4.res.rev,1,sum))
par(mfrow=c(1,1))


## choose burn-in

burn <- 10
its <- its-burn


## chop off the first burn iterations
ll.res.rev <- ll.res.rev[-(1:burn),]
ll2.res.rev <- ll2.res.rev[-(1:burn),]
ll3.res.rev <- ll3.res.rev[-(1:burn),]
ll4.res.rev <- ll4.res.rev[-(1:burn),]
res.rev$theta <- res.rev$theta[-(1:burn),]
res2.rev$theta <- res2.rev$theta[-(1:burn),]
res3.rev$theta <- res3.rev$theta[-(1:burn),]
res4.rev$theta <- res4.rev$theta[-(1:burn),]
par(mfrow=c(2,2))
plot.ts(apply(ll.res.rev,1,sum))
plot.ts(apply(ll2.res.rev,1,sum))
plot.ts(apply(ll3.res.rev,1,sum))
plot.ts(apply(ll4.res.rev,1,sum))
par(mfrow=c(1,1))


par(mfrow=c(3,2))
for(i in 1:K){
  plot(res.rev$theta[,i],col=1,type="l")
  lines(res2.rev$theta[,i],col=2)
  lines(res3.rev$theta[,i],col=3)
  lines(res4.rev$theta[,i],col=4)
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

pairs(res.rev$theta,pch=".")

## posterior correlations
crosscorr(as.mcmc(res.rev$theta))
crosscorr.plot(as.mcmc(res.rev$theta))


theta.rev.postmean <- apply(res.rev$theta,2,mean)
print(round(theta.rev.postmean,3))

## order the entities in terms of posterior mean (smallest to largest)
print(order(theta.rev.postmean))


## LOOIC

library(loo)

releff.rev <- relative_eff(exp(ll.res.rev),chain_id=rep(1,its))

loo.rev.res <- loo(ll.res.rev,r_eff=releff.rev)
print(loo.rev.res$estimates)
loo.rev.res$pointwise

## WAIC 
waic(ll.res.rev,r_eff=releff.rev)


#############################################
## EM algorithm


## set the seed
set.seed(2)

## maximum number of iterations
max.its <- 60

## stopping criteria for MSD - mean squared difference
stop.crit <- 10^(-16)

## initial values
theta.curr <- rbeta(K,a,b)

em.rev <- em.gpl.nb(K,y,s,w=ss$w,n.mat=n.mat,a,b,max.its,stop.crit,theta.curr)

theta.MAP.rev <- em.rev$theta[em.rev$its,]
print(round(theta.MAP.rev,3))

## order the entities in terms of MAP (smallest to largest)
print(order(theta.MAP.rev))

## posterior means and 95% equi-tailed credible intervals (plus MAP values)
pdf(file="pud-rev-postMAP.pdf",width=8,height=5,pointsize=14)
plot(1:K,theta.rev.postmean,ylim=c(0.3,0.6),pch=19,ylab=expression(theta[k]),xlab="k")
points(1:K,theta.MAP.rev,pch=4)
for(k in 1:K){
  lines(c(k,k),quantile(res.rev$theta[,k],c(0.025,0.975)))
}
dev.off()


####################################################################
##
## Davidson model
##
####################################################################


## read in the data 
puddings <- read.table("puddings.txt")
## columns are i j n_ij  w_ij  w_ji  t_ij

n.row <- dim(puddings)[1]


## convert data to y and s
n <- sum(puddings[,3])
K <- 6

y <- matrix(0,nrow=n,ncol=2)
s <- matrix(0,nrow=n,ncol=2)

k <- 0
for(m in 1:n.row){
  i <- puddings[m,1]
  j <- puddings[m,2]
  for(ell in 1:puddings[m,4]){
    k <- k+1
    y[k,1] <- i
    y[k,2] <- j
    s[k,1] <- 1
    s[k,2] <- 2
  }
  for(ell in 1:puddings[m,5]){
    k <- k+1
    y[k,1] <- j
    y[k,2] <- i
    s[k,1] <- 1
    s[k,2] <- 2
  }
  for(ell in 1:puddings[m,6]){
    k <- k+1
    y[k,1] <- i
    y[k,2] <- j
    s[k,1] <- 1
    s[k,2] <- 1
  }
}      

t <- s.to.t(s)



library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


source("gpl-stan-models.R")

## Davidson model

a <- 1
b <- 1
cc <- 1
d <- 1

Davidson_dat = list(K = K, n = n, y = y, t=t, a = a, b = b, c=cc, d=d)



set.seed(400)
fit_Davidson = stan(model_code = Davidson_code, data = Davidson_dat, iter = 5000, chains = 4)


##plot(fit_Davidson) 
##traceplot(fit_Davidson)
pairs(fit_Davidson) 
print(fit_Davidson,digits=3)
res.Dav = extract(fit_Davidson, permuted = TRUE)


dav.postmean.lambda <- apply(res.Dav$lambda,2,mean)
print(round(dav.postmean.lambda,3))
print(rev(order(round(dav.postmean.lambda,3))))

## posterior means and 95% equi-tailed credible intervals
pdf(file="pud-dav-post.pdf",width=8,height=5,pointsize=14)
plot(1:K,dav.postmean.lambda,ylim=c(0,3),pch=19,ylab=expression(lambda[k]),xlab="k")
for(k in 1:K){
  lines(c(k,k),quantile(res.Dav$lambda[,k],c(0.025,0.975)))
}
dev.off()

hist(res.Dav$delta,freq=FALSE)
print(mean(res.Dav$delta))
print(quantile(res.Dav$delta,c(0.025,0.975)))


library(loo)

ll.dav <- matrix(0,nrow=dim(res.Dav$lambda)[1],ncol=n)
for(i in 1:dim(res.Dav$lambda)[1]){  
  ll.dav[i,] <- Davidson.ll(y,t,res.Dav$lambda[i,],res.Dav$delta[i])
}

releff.dav <- relative_eff(exp(ll.dav),chain_id=rep(1,dim(res.Dav$lambda)[1]))


loo.dav <- loo(ll.dav,r_eff=releff.dav)
print(loo.dav)
waic.dav <- waic(ll.dav,r_eff=releff.dav)
waic.dav

