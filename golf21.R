## Analysis of 2021 PGA golf data from owgr.com

## author: daniel.henderson@newcastle.ac.uk

## read in GPL functions
source("gpl-functions.R")

## read in the results data
master <- read.table("golf21.txt")

## the players
players <- unique(master$Name)  ## factor, in order of appearance


## number of players
K <- length(players)

## possible outcomes
unique(master[,1])
## have to deal with WD = withdrawn, DQ = disqualified

## number of competitors in each tournament
ni <- as.numeric(table(master$id))
n <- length(ni)

## ordered list and ties set membership
##y <- matrix(0,nrow=n,ncol=K)  ## could reduce K to max(ni)
##s <- matrix(0,nrow=n,ncol=K)  ## could reduce K to max(ni)
y <- matrix(0,nrow=n,ncol=max(ni))  
s <- matrix(0,nrow=n,ncol=max(ni))
q <- ni
for(i in 1:n){
  y[i,1:ni[i]] <- match(master$Name[master$id==i],players)
  test <- unlist(strsplit(master[master$id==i,1],"T"))
  test <- as.numeric(test[test!=""])
  ranks <- rank(test,ties="min")
  s[i,1:ni[i]] <- match(ranks,unique(ranks))
  MC.WD.DQ <- which(master[master$id==i,1]=="MC" | master[master$i==i,1]=="WD" | master[master$i==i,1]=="DQ")
  if(length(MC.WD.DQ)>0){
    q[i] <- min(MC.WD.DQ)-1
    s[i,MC.WD.DQ] <- s[i,q[i]+1]
  }
}

image(y)
image(s)



## don't include the last tournament for now

y.all <- y
s.all <- s
q.all <- q
n.all <- length(q.all)

y <- y.all[-n.all,]
s <- s.all[-n.all,]
q <- q.all[-n.all]

n <- length(q)
K <- length(players)

t <- s.to.t(s)


## calculate summary statistics
pga.ss <- sum.stat.gpl(K,y,s,q)

## distribution of number of tournaments per golfer

pdf(file="uspga21-hist-tournaments.pdf",pointsize=14,width=8,height=5)
hist(pga.ss$cc,breaks=seq(0,35,1),xlab="Tournaments",main="")
dev.off()
barplot(table(pga.ss$cc))

table(pga.ss$cc)
print(median(pga.ss$cc))
print(sum(pga.ss$cc>=15))


## run Gibbs sampler

## beta prior hyperparameters
a <- rep(1,K)
b <- rep(1,K)


## number of iterations
its <- 10010


set.seed(1)
## initial values
theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
system.time(pga.res <- gs.gpl(K,y,s,v=pga.ss$v,w=pga.ss$w,delta=pga.ss$delta,a,b,its,theta.curr,z.curr))

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
pga.res2 <- gs.gpl(K,y,s,v=pga.ss$v,w=pga.ss$w,delta=pga.ss$delta,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
pga.res3 <- gs.gpl(K,y,s,v=pga.ss$v,w=pga.ss$w,delta=pga.ss$delta,a,b,its,theta.curr,z.curr)

theta.curr <- rbeta(K,a,b)
z.curr <- matrix(0,nrow=n,ncol=K) 
pga.res4 <- gs.gpl(K,y,s,v=pga.ss$v,w=pga.ss$w,delta=pga.ss$delta,a,b,its,theta.curr,z.curr)


## matrix of pointwise log-likelihoods over each datapoint
ll.res <- matrix(0,nrow=its,ncol=n)
ll2.res <- matrix(0,nrow=its,ncol=n)
ll3.res <- matrix(0,nrow=its,ncol=n)
ll4.res <- matrix(0,nrow=its,ncol=n)
for(i in 1:its){
  ll.res[i,] <- loglike.gpl(y,t,pga.res$theta[i,],m.i.star=pga.ss$m.i.star)
  ll2.res[i,] <- loglike.gpl(y,t,pga.res2$theta[i,],m.i.star=pga.ss$m.i.star)
  ll3.res[i,] <- loglike.gpl(y,t,pga.res3$theta[i,],m.i.star=pga.ss$m.i.star)
  ll4.res[i,] <- loglike.gpl(y,t,pga.res4$theta[i,],m.i.star=pga.ss$m.i.star)
}
ts.plot(ll.res,col=1:n)
ts.plot(ll2.res,col=1:n)
ts.plot(ll3.res,col=1:n)
ts.plot(ll4.res,col=1:n)

## decide on burn-in
burn <- 10
its <- its-burn

## chop off the first burn iterations
ll.res <- ll.res[-(1:burn),]
ll2.res <- ll2.res[-(1:burn),]
ll3.res <- ll3.res[-(1:burn),]
ll4.res <- ll4.res[-(1:burn),]
pga.res$theta <- pga.res$theta[-(1:burn),]
pga.res2$theta <- pga.res2$theta[-(1:burn),]
pga.res3$theta <- pga.res3$theta[-(1:burn),]
pga.res4$theta <- pga.res4$theta[-(1:burn),]
ts.plot(ll.res,col=1:n)
ts.plot(ll2.res,col=1:n)
ts.plot(ll3.res,col=1:n)
ts.plot(ll4.res,col=1:n)

par(mfrow=c(4,4))
for(i in 1:16){
  plot(pga.res$theta[,i],col=1,type="l")
  lines(pga.res2$theta[,i],col=2)
  lines(pga.res3$theta[,i],col=3)
  lines(pga.res4$theta[,i],col=4)
}
par(mfrow=c(1,1))


library(coda)

res.mat <- as.mcmc.list(list(as.mcmc(pga.res$theta),as.mcmc(pga.res2$theta),as.mcmc(pga.res3$theta),as.mcmc(pga.res4$theta)))

psrf <- gelman.diag(res.mat)
psrf$psrf
print(round(psrf$psrf,3))


## results from 1 chain

effectiveSize(pga.res$theta)
print(min(effectiveSize(pga.res$theta)))
which.min(effectiveSize(pga.res$theta))
acf(pga.res$theta[,1])
acf(pga.res$theta[,which.min(effectiveSize(pga.res$theta))])


pdf("uspga21-mix.pdf",width=15,height=5,pointsize=16)
par(mfrow=c(1,2))
plot(pga.res$theta[,which.min(effectiveSize(pga.res$theta))],type="l",xlab="Iteration",ylab=expression(theta[605]),main="")
acf(pga.res$theta[,which.min(effectiveSize(pga.res$theta))],main="")
par(mfrow=c(1,1))
dev.off()




## take into account those players who have played 15 or more tournaments

eligible.players <- which(pga.ss$cc>=15)


## top 10 
print(data.frame(1:10,players[eligible.players[rev(order(apply(pga.res$theta[,eligible.players],2,mean)))]][1:10],round(rev(sort(apply(pga.res$theta[,eligible.players],2,mean)))[1:10],3),round(apply(pga.res$theta[,eligible.players],2,quantile,0.025)[rev(order(apply(pga.res$theta[,eligible.players],2,mean)))[1:10]],3),round(apply(pga.res$theta[,eligible.players],2,quantile,0.975)[rev(order(apply(pga.res$theta[,eligible.players],2,mean)))[1:10]],3)))

## all eligible players in order of posterior mean
players[eligible.players[rev(order(apply(pga.res$theta[,eligible.players],2,mean)))]]



## pick out players in last tournament 
last.event <- y.all[n.all,1:q.all[n.all]]


post.mean.theta <- apply(pga.res$theta[,last.event],2,mean)
post.upper.theta <- apply(pga.res$theta[,last.event],2,quantile,0.975)
post.lower.theta <- apply(pga.res$theta[,last.event],2,quantile,0.025)


pdf(file="uspga21-hero-post.pdf",pointsize=12,width=8,height=5)
plot(rev(sort(post.mean.theta)),ylim=c(0,0.23),xlim=c(1,23),ylab=expression(theta),xlab="Model-based pre-tournament ranking",pch=19)
for(k in 1:length(post.mean.theta)){
  lines(c(k,k),c(post.lower.theta[rev(order(post.mean.theta))[k]],post.upper.theta[rev(order(post.mean.theta))[k]]),col=8)
  text(k,post.upper.theta[rev(order(post.mean.theta))[k]],players[last.event[rev(order(post.mean.theta))]][k],srt=45,pos=4,col=8)
}
dev.off()


## external odds

## https://www.cbssports.com/golf/news/2021-hero-world-challenge-odds-field-surprising-pga-picks-predictions-from-model-thats-nailed-7-majors/
odds <- c(33,10,12,9,12,18,16,20,18,25,16,11,16,20,8,14,22,28,22,66)

##Rory McIlroy 8-1
##Collin Morikawa 9-1
##Justin Thomas 10-1
##Viktor Hovland 11-1
##Xander Schauffele 12-1
##Bryson DeChambeau 12-1
##Jordan Spieth 14-1
##Sam Burns 16-1
##Tony Finau 16-1
##Scottie Scheffler 16-1
##Daniel Berger 18-1
##Webb Simpson 18-1
##Brooks Koepka 20-1
##Abraham Ancer 20-1
##Matthew Fitzpatrick 22-1
##Justin Rose 22-1
##Patrick Reed 25-1
##Tyrrell Hatton 28-1
##Harris English 33-1
##Henrik Stenson 66-1


## predictive simulations (with tie breaking for first)

sim.tournaments.y <- matrix(nrow=its,ncol=20)
sim.tournaments.s <- matrix(nrow=its,ncol=20)

winner <- rep(0,its)

for(i in 1:its){
  sim.t <- sim.data.gpl(1,pga.res$theta[i,last.event],t(as.matrix(1:20)))
  sim.tournaments.y[i,] <- last.event[sim.t$y]
  sim.tournaments.s[i,] <- sim.t$s
  tied.for.first <- sim.tournaments.y[i,sim.tournaments.s[i,]==1]
  n.tff <- length(tied.for.first)
  while(n.tff>1){
      print("tie breaker")
      ## play additional tie break holes until a single winner
      sim.t <- sim.data.gpl(1,pga.res$theta[i,tied.for.first],t(as.matrix(1:n.tff)))
      tied.for.first <- tied.for.first[sim.t$y][sim.t$s==1]
      n.tff <- length(tied.for.first) 
  }
  winner[i] <- tied.for.first
  print(c(i, winner[i]))
}

## final order
this.order <- c(18,15,14,12,10,4,16,5,2,17,9,11,20,3,1,6,19,13,7,8)

## for paper
odds.finalorder <- c(11,16,16,25,10,9,16,18,20,28,22,22,12,12,20,33,18,8,66,14)

## convert to probabilities and then normalise to sum to 1
round(1/(odds.finalorder+1)/sum(1/(odds.finalorder+1)),4)

## for table 5 in paper
print(cbind(players[last.event],cbind(table(players[winner])/its)[this.order],round(1/(odds.finalorder+1)/sum(1/(odds.finalorder+1)),4)))


## goodness of fit

## compare observed number of buckets (of ties) with simulated

## simulated
pdf(file="uspga21-hero-sets-ties.pdf",pointsize=12,width=8,height=5)
barplot(table(sim.tournaments.s[,20])/its,xlab="buckets",ylab="Probability")
dev.off()
## observed value
s.all[n.all,20]



