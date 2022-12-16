##################################################################
## Davidson model for paired comparisons

## author: daniel.henderson@newcastle.ac.uk

Davidson_code = '
  data {
    int<lower=1> K; // number of entities
    int<lower=1> n; // number of observations (paired comparisons)
    int<lower=0> y[n,2]; // paired comparisons
    int<lower=0,upper=1> t[n,1]; // ties indicator
    real<lower=0> a; // gamma prior parameter for lambda
    real<lower=0> b; // gamma prior parameter for lambda
    real<lower=0> c; // gamma prior parameter for delta
    real<lower=0> d; // gamma prior parameter for delta
 }
  parameters {
    real<lower=0> lambda[K]; 
    real<lower=0> delta;
  }
  model {
    real drll;
    lambda ~ gamma(a,b);
    delta ~ gamma(c,d);
    for(i in 1:n){
      drll = delta*sqrt(lambda[y[i,1]]*lambda[y[i,2]]);
      target += (1-t[i,1])*log(lambda[y[i,1]])+t[i,1]*log(drll)-log(lambda[y[i,1]]+lambda[y[i,2]]+drll); 
    }
  }
'
