
library(ggplot2)
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep.R")

tr = function(x) sum(diag(x))

stages = timepoints

genes = rownames(X)

timevector = rep(c("0h", "1h", "2h", "3h", "4h", "5h", 
                   "6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h"), each = 4)

stages_uniq = unique(timevector)

N = length(stages_uniq)

X = t(X)

dim(X)

p = ncol(X)

Y = vector("list", N)

j = 1

for(s in stages_uniq){
  
  ind = stringr::str_detect(stages, paste0("^", s))
  
  Y[[j]] = t(X[ind, ])
  
  j = j + 1
  
}

lapply(Y, dim)

# time varying scale-free glasso

pos = 1:N

S_t = kernel_cov_rep(X = Y, 
                     N = N,
                     pos = pos,
                     kernel = "epanechnikov")$S_t


nlambda = 50
lambda.min.ratio = 0.1

Theta_tvgl = array(0, c(p, p, N))

A_tvgl = vector("list", N)

for(i in 1:N){
  
  S_t_temp <- cov2cor(S_t[[i]])
  
  lambda.max = max(max(S_t_temp - diag(p)), 
                   -min(S_t_temp - diag(p)))
  lambda.min = lambda.min.ratio*lambda.max
  lambdas = exp(seq(log(lambda.max), 
                    log(lambda.min), 
                    length = nlambda))
  
  tvgl_res = list(icov = list(), 
                  df = c())
  
  j = 1
  
  for(lambda in lambdas){
    tvgl_res$icov[[j]] = glassoFast::glassoFast(S = S_t_temp,
                                                rho = lambda)$wi
    tvgl_res$df[j] = (sum(tvgl_res$icov[[j]] != 0) - ncol(tvgl_res$icov[[j]]))/2
    
    j = j + 1
    
  }
  
  # Select the optimal network using eBIC,
  
  loglik = rep(0, nlambda)
  for(j in 1:nlambda){
    loglik[j] = log(det(tvgl_res$icov[[j]])) - 
      tr(S_t_temp%*%tvgl_res$icov[[j]])
  }
  
  gamma = 0.5
  
  ebic = -N*loglik + log(N)*tvgl_res$df + 4*gamma*log(p)*tvgl_res$df

  plot(lambdas, ebic)
    
  ind = which.min(ebic)
  
  Theta_tvgl[,,i] = tvgl_res$icov[[ind]]
  
  A_temp = Theta_tvgl[,, i]
  
  A_temp = abs(A_temp)
  
  A_temp[A_temp != 0] = 1
  
  diag(A_temp) = 0
  
  A_tvgl[[i]] = A_temp
  
  cat("\r", i)
  
}

d_hat = d_hub = c()

Degree_all = matrix(0, N, p)

for(i in 1:N){
  
  G = graph_from_adjacency_matrix(A_tvgl[[i]],
                                  mode = "undirected",
                                  diag = FALSE)
  
  if(i == 1) init_h = which.max(degree(G))
  
  Degree_all[i, ] = degree(G)
  
  d_hat = c(d_hat, sum(degree(G)))
  
  d_hub = c(d_hub, degree(G)[init_h])
  
}

plot(pos, d_hat, type = "b")
abline(v = pos,
       lty = 2)

mtext(side = 3, 
      line = 0.5,
      at = pos,
      text = stages_uniq)

#

plot(pos, d_hub, type = "b")
abline(v = pos,
       lty = 2)

mtext(side = 3, 
      line = 0.5,
      at = pos,
      text = stages_uniq)


##

for(i in 1:N){
  
  A_temp = A_tvgl[[i]]
  
  ind = which(colSums(A_temp) == 0)
  
  A_temp = A_temp[-ind, -ind]
  
  G = graph_from_adjacency_matrix(A_temp,
                                  mode = "undirected",
                                  diag = FALSE)
  
  plot(G, 
       vertex.label = NA, 
       vertex.size = 10,
       main = stages_uniq[i])
  
}

