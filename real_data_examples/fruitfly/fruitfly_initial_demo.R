
library(mgm)
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep.R")
source("functions/kernel_cov.R")
source("functions/sf_glasso.R")

tr = function(x) sum(diag(x))

data("fruitfly_data")

X = fruitfly_data$data

stages = fruitfly_data$stages

genes = fruitfly_data$colnames

timevector = fruitfly_data$timevector

timevector

table(stages)

stages_uniq = unique(stages)

N = length(stages_uniq)
p = ncol(X)

Y = vector("list", N)

j = 1

for(s in stages_uniq){
  
  ind = which(stages == s)
  
  Y[[j]] = t(X[ind, ])
  
  j = j + 1
  
}

lapply(Y, dim)

# time varying scale-free glasso

pos = 1:4

S_t = kernel_cov_rep(X = Y, 
                     N = N,
                     pos = pos,
                     kernel = "gaussian")$S_t

nalpha = 50
alphas = seq(log(0.2), log(0.6), length.out = nalpha)
alphas = exp(alphas)

Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, N))

for(i in 1:N){
  
  # Select the optimal network using eBIC,
  
  Theta_temp = array(0, c(p, p, nalpha))
  
  gamma = 0.5
  
  loglik = ebic = rep(0, nalpha)
  
  j = 1
  
  for(alpha in alphas){
    
    Theta_temp[,, j] = sf_glasso(S = S_t[[i]], 
                                 alpha = alpha,
                                 K = 2) 
    
    loglik[j] = log(det(Theta_temp[,, j])) - 
      tr(S_t[[i]]%*%Theta_temp[,, j])
    
    df = (sum(Theta_temp[,, j] != 0) - 
            ncol(Theta_temp[,, j]))/2
    
    ebic[j] = -N*loglik[j] + log(N)*df + 4*gamma*log(p)*df
    
    j = j + 1
    
  }
  
  plot(alphas, ebic)
  
  ind = which.min(ebic)
  
  Theta_tvsfgl[,,i] = Theta_temp[,, ind]
  
  A_temp = Theta_tvsfgl[,, i]
  
  A_temp = abs(A_temp)
  
  A_temp[A_temp != 0] = 1
  
  diag(A_temp) = 0
  
  A_tvsfgl[,,i] = A_temp
  
  cat("\r", i)
  
}

d_hat = c()

for(i in 1:N){
  
  G = graph_from_adjacency_matrix(A_tvsfgl[,,i],
                                  mode = "undirected",
                                  diag = FALSE)
  
  d_hat = c(d_hat, sum(degree(G)))
  
  plot(G, 
       vertex.label = NA, 
       vertex.size = 10,
       main = stages_uniq[i])
  
}

plot(pos, d_hat, type = "b")
abline(v = pos,
       lty = 2)

mtext(side = 3, 
      line = 0.5,
      at = pos,
      text = stages_uniq)

# Smoothing from Zhou et al. (2008) when time-points
# are treated as equally spaced

N = nrow(X)
pos = 1:N

dim(X)

S_t = kernel_cov(X = t(X), 
                 pos = pos,
                 kernel = "gaussian")$S_t

nalpha = 50
alphas = seq(log(0.1), log(0.5), length.out = nalpha)
alphas = exp(alphas)

Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, N))

for(i in 1:N){
  
  # Select the optimal network using eBIC,
  
  Theta_temp = array(0, c(p, p, nalpha))
  
  gamma = 0.5
  
  loglik = ebic = rep(0, nalpha)
  
  j = 1
  
  for(alpha in alphas){
    
    Theta_temp[,, j] = sf_glasso(S = S_t[,,i], 
                                 alpha = alpha,
                                 K = 2) 
    
    loglik[j] = log(det(Theta_temp[,, j])) - 
      tr(S_t[,,i]%*%Theta_temp[,, j])
    
    df = (sum(Theta_temp[,, j] != 0) - 
            ncol(Theta_temp[,, j]))/2
    
    ebic[j] = -N*loglik[j] + log(N)*df + 4*gamma*log(p)*df
    
    j = j + 1
    
  }
  
  ind = which.min(ebic)
  
  Theta_tvsfgl[,,i] = Theta_temp[,, ind]
  
  A_temp = Theta_tvsfgl[,, i]
  
  A_temp = abs(A_temp)
  
  A_temp[A_temp != 0] = 1
  
  diag(A_temp) = 0
  
  A_tvsfgl[,,i] = A_temp
  
  cat("\r", i)
  
}

d_hat = c()

for(i in 1:N){
  
  G = graph_from_adjacency_matrix(A_tvsfgl[,,i],
                                  mode = "undirected",
                                  diag = FALSE)
  
  d_hat = c(d_hat, sum(degree(G)))
  
  plot(G, 
       vertex.label = NA, 
       vertex.size = 10,
       main = paste0(stages[i], 
                     "_time = ", timevector[i]))
  
}

stages_end = c()

j = 1

for(i in stages_uniq){
  
  stages_end[j] = tail(which(stages %in% i), n = 1)
  
  j = j + 1
  
}

plot(1:N, d_hat, type = "b")
abline(v = stages_end, 
       lty = 2, 
       lwd = 2, 
       col = "red")

mtext(side = 3, 
      line = 0.5,
      at = stages_end,
      text = stages_uniq)
