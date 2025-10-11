
library(mgm)
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/sf_glasso.R")
source("real_data_examples/fruitfly/weighted_cov.R")
source("real_data_examples/fruitfly/weighted_cov_timeseries.R")

tr = function(x) sum(diag(x))

data("fruitfly_data")

X = fruitfly_data$data

stages = fruitfly_data$stages

genes = fruitfly_data$colnames

timevector = fruitfly_data$timevector

timevector

vmin = min(timevector)
vmax = max(timevector)

s = (timevector - vmin)/vmax

N = length(timevector)

plot(s, 
     y = rep(0, N), 
     ylim = c(0, 1),
     pch = " ")

abline(v = s, lty = 2)

table(stages)

stages_uniq = unique(stages)

stages_end = c()

j = 1

for(i in stages_uniq){
  
  stages_end[j] = tail(which(stages %in% i), n = 1)
  
  j = j + 1
  
}

stages_end

abline(v = s[stages_end], 
       lty = 2, 
       col = "red",
       lwd = 2)

p = ncol(X)

Y = vector("list", length = length(stages_uniq))

j = 1

for(u in stages_uniq){
  
  ind = which(stages == u)
  
  Y[[j]] = t(X[ind, ])
  
  j = j + 1
  
}

lapply(Y, dim)

# time varying scale-free glasso

# First, aggregate values from each phase of the
# fruitfly life cycle

pos = s[stages_end]

N = length(Y)

#h <- 5.848/nrow(X)^(1/3)

S_t = weighted_cov(Y = Y, 
                   s = pos)

lapply(S_t, dim)

n = lapply(Y, function(x) dim(x)[2])

n = unlist(n)

nalpha = 40
alphas = seq(log(0.3), log(0.4), length.out = nalpha)
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

d_hat

plot(pos, d_hat, type = "b")
abline(v = s[stages_end], 
       lty = 2, 
       lwd = 2, 
       col = "red")

# fixed alpha

alpha = 0.21

Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, N))

for(i in 1:N){
  
  Theta_temp = sf_glasso(S = S_t[[i]], 
                         alpha = alpha,
                         K = 4) 
  
  Theta_tvsfgl[,,i] = Theta_temp
  
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

d_hat

plot(pos, d_hat, type = "b")
abline(v = s[stages_end], 
       lty = 2, 
       lwd = 2, 
       col = "red")

mtext(side = 3, 
      line = 0.5,
      at = s[stages_end],
      text = stages_uniq)

# Zhou et al (2008) smoothing

dim(X)

pos = seq(0, 1, length.out = 20)

S_t = weighted_cov_timeseries(X = t(X), 
                              s = s,
                              pos = pos,
                              kernel = "gaussian")$S_t

dim(S_t)

N = nrow(X)

nalpha = 30
alphas = seq(log(0.2), log(0.35), length.out = nalpha)
alphas = exp(alphas)

Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, length(pos)))

for(i in 1:length(pos)){
  
  # Select the optimal network using eBIC,
  
  # Note! This model selection criterion is probably
  # not suitable in the case when there are no 
  # repeated samples at each time step
  
  Theta_temp = array(0, c(p, p, nalpha))
  
  gamma = 0.5
  
  loglik = ebic = rep(0, nalpha)
  
  j = 1
  
  for(alpha in alphas){
    
    Theta_temp[,, j] = sf_glasso(S = S_t[,,i], 
                                 alpha = alpha,
                                 K = 4) 
    
    loglik[j] = log(det(Theta_temp[,, j])) - 
      tr(S_t[,,i]%*%Theta_temp[,, j])
    
    df = (sum(Theta_temp[,, j] != 0) - 
            ncol(Theta_temp[,, j]))/2
    
    ebic[j] = -N*loglik[j] + 
      log(N)*df + 4*gamma*log(p)*df
    
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

for(i in 1:length(pos)){
  
  G = graph_from_adjacency_matrix(A_tvsfgl[,,i],
                                  mode = "undirected",
                                  diag = FALSE)
  
  d_hat = c(d_hat, sum(degree(G)))
  
  plot(G, 
       vertex.label = NA, 
       vertex.size = 10,
       main = paste0(stages[i], 
                     "time=", round(pos[i], 3)))
  
}

plot(pos, d_hat, type = "l")
abline(v = pos, 
       lty = 2, 
       lwd = 1)

abline(v = s[stages_end],
       lty = 2,
       col = "red",
       lwd = 2)

# Fixed lambda value,

alpha = alphas[ind]

Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, length(pos)))

for(i in 1:length(pos)){
  
  # fixed alpha value,
  
  Theta_temp = sf_glasso(S = S_t[,,i], 
                                 alpha = alpha,
                                 K = 4) 
  
  Theta_tvsfgl[,,i] = Theta_temp
  
  A_temp = Theta_tvsfgl[,, i]
  
  A_temp = abs(A_temp)
  
  A_temp[A_temp != 0] = 1
  
  diag(A_temp) = 0
  
  A_tvsfgl[,,i] = A_temp
  
  cat("\r", i)
  
}

d_hat = c()

for(i in 1:length(pos)){
  
  G = graph_from_adjacency_matrix(A_tvsfgl[,,i],
                                  mode = "undirected",
                                  diag = FALSE)
  
  d_hat = c(d_hat, sum(degree(G)))
  
  plot(G, 
       vertex.label = NA, 
       vertex.size = 10,
       main = paste0(stages[i], 
                     "time=", round(pos[i], 3)))
  
}

plot(pos, d_hat, type = "l")
abline(v = pos, 
       lty = 2, 
       lwd = 1)

abline(v = s[stages_end], 
       lty = 2, 
       lwd = 2, 
       col = "red")

mtext(side = 3, 
      line = 0.5,
      at = s[stages_end],
      text = stages_uniq)

##

# Aggregate, say, 10 observed values
# of the fruitfly life cycle as one dynamic
# covariance

source("functions/kernel_cov_rep.R")

m = matrix(c(1:67,rep(NA, 3)), 
           nrow = 10, 
           ncol = 7, 
           byrow = F)

Y = vector("list", length = ncol(m))

for(j in 1:ncol(m)){
  
  x = m[, j]
  
  x = x[!is.na(x)]
  
  Y[[j]] = t(X[x, ])
  
}

lapply(Y, dim)

N = length(Y)

S_t = kernel_cov_rep(X = Y,
                     N = N,
                     pos = 1:N)$S_t

lapply(S_t, dim)

n = lapply(Y, function(x) dim(x)[2])

n = unlist(n)

nalpha = 40
alphas = seq(log(0.2), log(0.33), length.out = nalpha)
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
  
  plot(alphas,
       ebic,
       main = paste0("eBIC at ", i))
  
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
       main = paste0(i))
  
}

d_hat

d = seq(10, 60, by = 10)

d = c(d, 67)

plot(s[d], d_hat, type = "l")
abline(v = s[d], 
       lty = 2, 
       lwd = 1)

abline(v = s[stages_end], 
       lty = 2, 
       lwd = 2, 
       col = "red")

mtext(side = 3, 
      line = 0.5,
      at = s[stages_end],
      text = stages_uniq)

# Fixed alpha

alpha = 0.2

Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, N))

for(i in 1:N){
    
  Theta_tvsfgl[,, i] = sf_glasso(S = S_t[[i]], 
                                 alpha = alpha,
                                 K = 2) 
  
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
       main = paste0(i))
  
}

d_hat

d = seq(10, 60, by = 10)

d = c(d, 67)

plot(s[d], d_hat, type = "l")
abline(v = s[d], 
       lty = 2, 
       lwd = 1)

abline(v = s[stages_end], 
       lty = 2, 
       lwd = 2, 
       col = "red")

mtext(side = 3, 
      line = 0.5,
      at = s[stages_end],
      text = stages_uniq)
