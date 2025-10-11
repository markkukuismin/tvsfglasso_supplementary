
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep.R")
source("functions/sf_glasso.R")
source("functions/binary_class.R")
source("functions/loggle.graph.R")

tr = function(x) sum(diag(x))

# Simulate time-varying graph data using 
# loggle.graph from Yang & Peng (2018)

set.seed(1)

p = 100
n = 10 # Nmb of repetitions, 10, 50, 100
N = 200 # nmb of time points

# Repeat the process SM times to take into account
# the random error caused by the simulation process

SM = 10

Results = vector("list", SM)

example <- loggle.graph(p = p, n = n, N = N)

length(example$Omega.true)

Theta_true = example$Omega.true

A_true = lapply(Theta_true, as.matrix)

Sigma_true = lapply(Theta_true, as.matrix)

for(i in 1:N){
  
  A_temp = A_true[[i]]
  
  A_temp[A_temp != 0] = 1
  
  A_true[[i]] = A_temp
  
  Sigma_true[[i]] = Cholesky(Theta_true[[i]])
  
}

# Use different bandwidth parameters,

hband = c(0.01, 0.5, 1, 1.5, 2.5)

for(h in hband){
  
  for(sim_i in 1:SM){
    
    Y <- vector("list", N)
    
    for(i in 1:N){
      
      Y[[i]] <- t(rmvn.sparse(n = n, mu = rep(0, p), 
                              Sigma_true[[i]]))
      
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
    
    Frob_norm = lapply(S_t, function(x) norm(x - S, type = "F"))
    
    Frob_norm = unlist(Frob_norm)
    
    A_tvgl = array(0, c(p, p, N))
    
    R_tvgl = matrix(0, 12, 1)
    
    lambda.min.ratio = 0.1
    
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
      
      df = tvgl_res$df
      gamma = 0.5
      
      ebic = -N*loglik + log(N)*df + 4*gamma*log(p)*df
      
      ind = which.min(ebic)
      
      #A_temp = tvgl_res$path[[ind]]
      A_temp = tvgl_res$icov[[ind]]
      
      A_temp[A_temp != 0] = 1
      
      diag(A_temp) = 0
      
      R_tvgl = cbind(R_tvgl,
                     unlist(binary_class(A = A_true[[i]], 
                                         A_est = A_temp))
      )
      
      A_tvgl[,,i] = A_temp
      
    }
    
    R_tvgl = R_tvgl[, -1]
    
    # time-varying scale-free glasso
    
    Theta_tvsfgl = A_tvsfgl = array(0, c(p, p, N))
    
    R_tvsfgl = matrix(0, 12, 1)
    
    nalpha = 25
    alphas = seq(log(0.05), log(0.3), length.out = nalpha)
    alphas = exp(alphas)
    
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
    
    m = rep(c("tvglasso", "tvsfglasso"), each = N)
    
    D = data.frame(TP = c(R_tvgl["TP", ],
                          R_tvsfgl["TP", ]),
                   FP = c(R_tvgl["FP", ],
                          R_tvsfgl["FP", ]),
                   FN = c(R_tvgl["FN", ],
                          R_tvsfgl["FN", ]),
                   TN = c(R_tvgl["TN", ],
                          R_tvsfgl["TN", ]),
                   Pre = c(R_tvgl["Precision", ],
                           R_tvsfgl["Precision", ]),
                   TPR = c(R_tvgl["TPR", ],
                           R_tvsfgl["TPR", ]),
                   MCC = c(R_tvgl["MCC", ],
                           R_tvsfgl["MCC", ]),
                   FDR = c(R_tvgl["FDR", ],
                           R_tvsfgl["FDR", ]),
                   FPR = c(R_tvgl["FPR", ],
                           R_tvsfgl["FPR", ]),
                   F1 = c(R_tvgl["F1", ],
                          R_tvsfgl["F1", ]),
                   JI = c(R_tvgl["JI", ],
                          R_tvsfgl["JI", ]),
                   ED = c(R_tvgl["ED", ],
                          R_tvsfgl["ED", ]),
                   Method = m,
                   Sim_round = sim_i)
    
    
    Results[[sim_i]] = D
    
    #cat("\r", sim_i)
    
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
  
  fn = paste0("simulations/simulations_with_more_rounds/results/yang_peng/tvglasso_and_tvsfglasso_bandwidth_res/yang_peng_mod_", n, "_", p, "_", N, "_", h, ".txt")
  write.table(Sim_res, file = fn, row.names = FALSE)
  
  cat("\r", h)
  
}