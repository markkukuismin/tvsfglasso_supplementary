
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep.R")
source("functions/sf_glasso.R")
source("functions/binary_class.R")
source("functions/generate_tv_sf_data.R")

tr = function(x) sum(diag(x))

# Simulate time-varying graph data using 
# huge_generate_tv_sf_data

set.seed(1)

p = 100
n = 10 # Nmb of repetitions, 10, 50, 100
N = 200 # nmb of time points

# Repeat the process SM times to take into account
# the random error caused by the simulation process

SM = 10

Results = vector("list", SM)

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

# Use different bandwidth parameters,

h = 1

ebic_gamma = c(0.1, 0.2, 0.5, 0.7, 0.9)

for(gamma in ebic_gamma){
  
  for(sim_i in 1:SM){
    
    Y = vector("list", N)
    
    for(i in 1:N){
      
      Y[[i]] = t(MASS::mvrnorm(n, mu = rep(0, p), Sigma = example$Sigma.true[[i]]))
      
    }
    
    nlambda = 25
    
    # time-varying glasso
    
    pos = 1:N
    
    S_t = kernel_cov_rep(X = Y, 
                         N = N,
                         pos = pos,
                         h = h,
                         kernel = "epanechnikov")$S_t
    
    S = cov(t(Reduce(cbind, Y)))
    
    lambda.min.ratio = 0.1
    
    # time-varying scale-free glasso
    
    Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, N))
    
    R_tvsfgl = matrix(0, 12, 1)
    
    nalpha = 25
    alphas = seq(log(0.05), log(0.3), length.out = nalpha)
    alphas = exp(alphas)
    
    for(i in 1:N){
      
      # Select the optimal network using eBIC,
      
      Theta_temp = array(0, c(p, p, nalpha))
      
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
      
      ind = which.min(ebic)
      
      Theta_tvsfgl[,,i] = Theta_temp[,, ind]
      
      A_temp = Theta_tvsfgl[,, i]
      
      A_temp = abs(A_temp)
      
      A_temp[A_temp != 0] = 1
      
      diag(A_temp) = 0
      
      R_tvsfgl = cbind(R_tvsfgl,
                       unlist(binary_class(A = A_true[[i]], 
                                           A_est = A_temp))
      )
      
      A_tvsfgl[,,i] = A_temp
      
    }
    
    R_tvsfgl = R_tvsfgl[, -1]
    
    m = rep("tvsfglasso", each = N)
    
    D = data.frame(TP = R_tvsfgl["TP", ],
                   FP = R_tvsfgl["FP", ],
                   FN = R_tvsfgl["FN", ],
                   TN = R_tvsfgl["TN", ],
                   Pre = R_tvsfgl["Precision", ],
                   TPR = R_tvsfgl["TPR", ],
                   MCC = R_tvsfgl["MCC", ],
                   FDR = R_tvsfgl["FDR", ],
                   FPR = R_tvsfgl["FPR", ],
                   F1 = R_tvsfgl["F1", ],
                   JI = R_tvsfgl["JI", ],
                   ED = R_tvsfgl["ED", ],
                   Method = m,
                   Sim_round = sim_i)
    
    
    Results[[sim_i]] = D
    
  }
  
  Sim_res = dplyr::bind_rows(Results)
  
  Sim_res = tidyr::pivot_longer(Sim_res, 
                                cols = 1:12, 
                                values_to = "Value",
                                names_to = "Metric")
  
  Sim_res$N = N
  Sim_res$n = n
  Sim_res$p = p
  Sim_res$h = h
  
  ##
  
  fn = paste0("simulations/simulations_with_more_rounds/ebic_sensitivity/liu_ihler_mod_", n, "_", p, "_", N, "_", gamma, ".txt")
  write.table(Sim_res, file = fn, row.names = FALSE)
  
  cat("\r", gamma)
  
}