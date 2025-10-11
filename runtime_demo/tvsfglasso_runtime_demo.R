
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep.R")
source("functions/sf_glasso.R")
source("functions/huge_generate_tv_sf_data.R")

tr = function(x) sum(diag(x))

set.seed(1)

p = 500 # nmb of variables
n = 10 # nmb of repetitions
N = 200 # nmb of time points

example <- huge_generate_tv_sf_data(p = p, 
                                    n = n, 
                                    N = N, 
                                    alpha = 0.2)

Y = example$X

Y = vector("list", N)
  
for(i in 1:N){
    
  Y[[i]] = t(MASS::mvrnorm(n, 
                           mu = rep(0, p), 
                           Sigma = example$Sigma.true[[i]]))
    
}
  
nalpha = 25
alphas = seq(log(0.05), log(0.3), length.out = nalpha)
alphas = exp(alphas)

pos = 1:N

start_time = proc.time()
  
S_t = kernel_cov_rep(X = Y, 
                     N = N,
                     pos = pos,
                     kernel = "epanechnikov")$S_t
    
for(i in 1:N){
  
  for(alpha in alphas){
    
    Theta_temp = sf_glasso(S = S_t[[i]], 
                           alpha = alpha,
                           K = 2) 
    
  }
  
}

stop_time = proc.time()

stop_time - start_time

#   user  system elapsed 
# 123.07   36.05  168.34 