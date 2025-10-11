
source("functions/TVGL_hallac.R")

library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/binary_class.R")
source("functions/generate_tv_sf_data.R")

tr = function(x) sum(diag(x))

# Simulate time-varying graph data using 
# huge_generate_tv_sf_data

set.seed(1)

p = 100
n = 50 # Nmb of repetitions, 10, 50, 100
N = 100 # nmb of time points

# Repeat the process SM times to take into account
# the random error caused by the simulation process

# SM = 10
# 
# Results = vector("list", SM)

example <- generate_tv_sf_data(p = p, 
                               n = n, 
                               N = N, 
                               alpha = 0.2)

Y = example$X

Theta_true = example$Omega.true

A_true = lapply(Theta_true, as.matrix)

for(i in 1:N){
  
  A_temp = A_true[[i]]
  
  A_temp[A_temp != 0] = 1
  
  A_true[[i]] = A_temp
  
}


Y <- vector("list", N)
  
for(i in 1:N){
    
  Y[[i]] <- t(MASS::mvrnorm(n, mu = rep(0, p), Sigma = example$Sigma.true[[i]]))
    
}

fun = function(x) t(x)

Y <- lapply(Y, fun)

length(Y)
dim(Y[[1]])

system.time(res <- TVGL(Y, lambda1 = 0.5, lambda2 = 0.5))

