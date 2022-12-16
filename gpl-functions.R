##
## R functions for results in 
##
## "Modelling and analysis of rank ordered data with ties via a
## generalized Plackett-Luce model"
##
## author: daniel.henderson@newcastle.ac.uk

rg <- function(n,p)
{
  ## sample n random variates from the Geom(p) distribution on {1,2,...}

  1+rgeom(n,p)
}

loglike.gpl <- function(y,t,theta,m.i.star)
{
  ## Generalized (geometric) Plackett-Luce (GPL)
  ## loglikelihoods over multiple datasets

  n = max(dim(y)[1],1)
  if(n == 1){
    y = t(matrix(y))
    t = t(matrix(t))
  }
  ll.vec = rep(0,n)
  for(i in 1:n){
    n.i = length(y[i,y[i,]>0])
    for(j in 1:(m.i.star[i])){
      ll.vec[i] = ll.vec[i]+log(theta[y[i,j]])-log(1-exp(sum(log(1-theta[y[i,j:n.i]]))))+(1-t[i,j])*sum(log(1-theta[y[i,(j+1):n.i]])) + t[i,j]*log(1-exp(sum(log(1-theta[y[i,(j+1):n.i]]))))
    }
  }
  return(ll.vec)
}


sim.data.gpl <- function(n,theta,these.entities,q=rep(length(theta),n))
{
  ## simulate data from the GPL model 
  ## simulate latent variables, compute ranks, assign set membership
  ## compute number of ordered sets

  K <- length(theta)                ## number of entities
  U <- matrix(0,nrow=n,ncol=K)      ## latent variables
  r.U <- matrix(0,nrow=n,ncol=K)    ## ranks
  SM <- matrix(0,nrow=n,ncol=K)     ## set membership
  s <- matrix(0,nrow=n,ncol=K)      ## sorted set membership
  y <- matrix(0,nrow=n,ncol=K)      ## ordering of entities
  n.i <- rep(0,n)                   ## number of entities in ith comparison
  m.i <- rep(K,n)                   ## top m_i
  m.i.star <- rep(0,n)              ## min(m.i,n.i-1)
  v <- rep(0,n)                     ## number of ordered sets containing 
                                    ## the top m_i entities but not including 
                                    ## the s_in_ith ordered set if it does not 
                                    ## involve a tie
  for(i in 1:n){
    n.i[i] <- sum(these.entities[i,]>0)
    U[i,1:n.i[i]] <- rg(n.i[i],theta[these.entities[i,these.entities[i,]>0]])
    r.U[i,1:n.i[i]] <- rank(U[i,1:n.i[i]])
    SM[i,1:n.i[i]] <- match(r.U[i,1:n.i[i]],sort(unique(r.U[i,1:n.i[i]])))
    s[i,1:n.i[i]] <- sort(SM[i,1:n.i[i]])
    y[i,1:n.i[i]] <- these.entities[i,1:n.i[i]][order(SM[i,1:n.i[i]])]
  }

  ## post-processing for top q rankings
  s.true <- s
  for(i in 1:n){
    q[i] <- min(q[i],n.i[i])
    m.i[i] <- sum(s.true[i,1:n.i[i]]<= s.true[i,q[i]])
    m.i.star[i] <- min(m.i[i],n.i[i]-1)
    v[i] <- s.true[i,m.i.star[i]]
    if(m.i[i]<n.i[i]){
      s[i,(m.i[i]+1):n.i[i]] <- s.true[i,m.i[i]]+1
    }
  }
  return(list(y=y,s=s))
}


sum.stat.gpl <- function(K,y,s,q=rep(K,dim(y)[1]))
{

  ## calculate summary statistics based on y and s (and q)
  
  n <- dim(y)[1]
  n.i <- apply(y>0,1,sum)           ## number of entities in each comparison
  m.i <- rep(K,n)                   ## top m_i
  m.i.star <- rep(0,n)              ## min(m.i,n.i-1)
  nos <- rep(0,n)                   ## number of ordered sets
  v <- rep(0,n)                     ## number of ordered sets containing 
                                    ## the top m_i entities but not including 
                                    ## the s_in_ith ordered set if it does not 
                                    ## involve a tie

  for(i in 1:n){
    q[i] <- min(q[i],n.i[i])
    m.i[i] <- sum(s[i,1:n.i[i]]<= s[i,q[i]])
    m.i.star[i] <- min(m.i[i],n.i[i]-1)
    nos[i] <- max(s[i,1:n.i[i]])
    v[i] <- s[i,m.i.star[i]]
  }

  ## number of comparisons involving each entity
  cc <- rep(0,K)

  for(k in 1:K){
    cc[k] <- sum(y==k)
  }

  ## compute delta_{ijk}
  delta <- array(0,c(n,K,K))
  for(i in 1:n){
    for(j in 1:s[i,n.i[i]]){
      for(k in 1:K){
        delta[i,j,k] <- sum(y[i,s[i,1:n.i[i]]>=j]==k)
      }
    }
  }
  
  ## number of sets ranked higher than the set containing entity k
  d <- rep(0,K)
  for(k in 1:K){
    for(i in 1:n){
      if(s[i,n.i[i]]>1){
        d[k] <- d[k]+sum(delta[i,2:s[i,n.i[i]],k])
      }
    }
  }

  ## number of comparisons in which entity k appears in the top m.i
  ## (and is not in position n.i on its own).
  w <- rep(0,K)
  for(k in 1:K){
    for(i in 1:n){
      w[k] <- w[k]+sum(delta[i,1:v[i],k])
    }
    w[k] <- w[k]-d[k]
  }

  return(list(delta=delta,v=v,w=w,d=d,n.i=n.i,cc=cc,m.i=m.i,m.i.star=m.i.star,q=q))
}

loglike.ss.gpl <- function(ss,theta)
{
  ## GPL model 
  ## loglikelihood based on sufficient statistics
  
  ll <- sum(ss$w*log(theta)+ss$d*log(1-theta))
  n <-  length(ss$n.i)
  for(i in 1:n){
    for(j in 1:ss$v[i]){ 
      ll <- ll-log(1-exp(sum(log(1-theta[ss$delta[i,j,]>0]))))
    }
  }
  return(ll)  
}


gs.gpl <- function(K,y,s,v,w,delta,a,b,its,theta.curr,z.curr)
{
   ## Gibbs sampler for GPL model
   
   n <- dim(y)[1]

   ## store samples
   theta.sim <- matrix(0,nrow=its,ncol=K)
   z.sim <- array(0,c(its,n,K))

   for(ell in 1:its){
   ##print(ell)

     for(i in 1:n){
       for(j in 1:v[i]){
         ##print(c(i,j,1-exp(sum(log(1-theta.curr[y[i,s[i,]>=j]])))))
         z.curr[i,j] <- rg(1,1-exp(sum(log(1-theta.curr[y[i,s[i,]>=j]]))))
       }
     }

     ## compute zeta.curr
     zeta.curr <- rep(0,K)
     for(k in 1:K){
       zeta.curr[k] <- sum(delta[,,k]*z.curr)
     }

     for(k in 1:K){
       ##print(c(ell,k,b[k]+zeta.curr[k]-w[k]))
       theta.curr[k] <- rbeta(1,a[k]+w[k],b[k]+zeta.curr[k]-w[k])
     }

     theta.sim[ell,] <- theta.curr
     z.sim[ell,,] <- z.curr
   }
   return(list(theta=theta.sim,z=z.sim))
}

em.gpl <- function(K,y,s,v,w,delta,a,b,max.its,stop.crit=10^(-16),theta.curr)
{
   ## EM algorithm for GPL model
   
   n <- dim(y)[1]

   ## store estimates
   theta.est <- matrix(0,nrow=max.its,ncol=K)

   obj <- 10^6
   ell <- 0
   
   while(ell < max.its & obj > stop.crit){
     ell <- ell+1

     theta.old <- theta.curr

     ## compute rho.curr
     rho.curr <- rep(0,K)
     for(k in 1:K){
       for(i in 1:n){
         for(j in 1:v[i]){

           rho.curr[k] <- rho.curr[k]+(delta[i,j,k]/(1-exp(sum(log(1-theta.curr[y[i,s[i,]>=j]])))))
         }
       }       
     }

     for(k in 1:K){
       theta.curr[k] <- (a[k]+w[k]-1)/(a[k]+b[k]+rho.curr[k]-2)
     }

     obj <- mean((theta.old-theta.curr)^2)
     ##print(obj)

     theta.est[ell,] <- theta.curr
   }
   return(list(theta=theta.est,its=ell,obj=obj))
}



t.to.s <- function(t)
{
  ## GPL model
  ## determine s from t
  
  K <- dim(t)[2]+1
  n <- dim(t)[1]
  s <- matrix(0,nrow=n,ncol=K)
  for(i in 1:n){
    s[i,] <- 1+((1:K)-cumsum(c(1,t[i,])))
  }
  s
}

s.to.t <- function(s)
{

  ## GPL model
  ## determine t from s
  
  K <- dim(s)[2]
  n <- dim(s)[1]
  t <- matrix(0,nrow=n,ncol=K-1)
  for(i in 1:n){
    t[i,] <- ifelse(diff(s[i,])==0,1,0)
  }
  t
}


ltp <- function(theta.x,theta.y)
{
  ## GPL model
  ## probability that W_x < W_y when W_i ~ Geom(theta_i) indep.
  
  theta.x*(1-theta.y)/(1-(1-theta.x)*(1-theta.y))
}

eqp <- function(theta.x,theta.y)
{
  ## GPL model
  ## probability that W_x = W_y when W_i ~ Geom(theta_i) indep.
  
  theta.x*theta.y/(1-(1-theta.x)*(1-theta.y))
}




######################################################################

## alternative negative binomial latent variable formulation

rnb <- function(n,r,p)
{
  ## sample n random variates from the NBinom(r,p) distribution on
  ## {r,r+1,...}, for r=1,2,...

  replicate(n,sum(rg(r,p)))	
}


gs.gpl.nb <- function(K,y,s,w,n.mat,a,b,its,theta.curr,z.curr)
{
   ## Gibbs sampler for paired comparisons
   ## Using alternative negative binomial latent variable formulation
   
   n <- dim(y)

   ## store samples
   theta.sim <- matrix(0,nrow=its,ncol=K)
   z.sim <- array(0,c(its,K,K))

   for(ell in 1:its){
     for(i in 1:(K-1)){
       for(j in (i+1):K){
         if(n.mat[i,j]>0){
           z.curr[i,j] <- rnb(1,n.mat[i,j],1-(1-theta.curr[i])*(1-theta.curr[j]))
         }
       }
     }

     xi.curr <- rep(0,K)
     for(i in 1:K){
       xi.curr[i] <- sum(z.curr[,i])+sum(z.curr[i,])
     }

     for(i in 1:K){
       theta.curr[i] <- rbeta(1,a[i]+w[i],b[i]+xi.curr[i]-w[i])
     }

     theta.sim[ell,] <- theta.curr
     z.sim[ell,,] <- z.curr
   }
   return(list(theta=theta.sim,z=z.sim))
}

em.gpl.nb <- function(K,y,s,w,n.mat,a,b,max.its,stop.crit=10^(-16),theta.curr)
{
   ## EM algorithm for paired comparisons
   ## Using alternative negative binomial latent variable formulation
   
   n <- dim(y)

   ## store samples
   theta.est <- matrix(0,nrow=max.its,ncol=K)

   obj <- 10^6
   ell <- 0
   

   while(ell < max.its & obj > stop.crit){
     ell <- ell+1

     theta.old <- theta.curr

     rho.curr <- rep(0,K)
     
     for(i in 1:K){
       ## OK for sums to include i=j as n_{ii}=0
       rho.curr[i] <- sum(n.mat[i,i:K]/(1-(1-theta.curr[i])*(1-theta.curr[i:K])))+sum(n.mat[1:i,i]/(1-(1-theta.curr[1:i])*(1-theta.curr[i])))
     }     

 
     for(i in 1:K){
       theta.curr[i] <- (a[i]+w[i]-1)/(a[i]+b[i]+rho.curr[i]-2)
     }

     obj <- mean((theta.old-theta.curr)^2)
     ##print(obj)

     theta.est[ell,] <- theta.curr
   }
   return(list(theta=theta.est,its=ell,obj=obj))
}



Davidson.ll <- function(y,t,lambda,delta)
{
  ## calculate log likelihood under Davidson (1970) model
  n <- dim(y)[1]
  ll <- rep(0,n)
  for(i in 1:n){
    drll = delta*sqrt(lambda[y[i,1]]*lambda[y[i,2]])
    ll[i] <- (1-t[i,1])*log(lambda[y[i,1]])+t[i,1]*log(drll)-log(lambda[y[i,1]]+lambda[y[i,2]]+drll)
  }
  ll
}


