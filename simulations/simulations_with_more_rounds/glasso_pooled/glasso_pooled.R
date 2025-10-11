
library(igraph)
library(huge)
library(hglasso)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep.R")
source("functions/sf_glasso.R")
source("functions/binary_class.R")
source("functions/hglasso_generate_tv_sf_data.R")
source("functions/huge_generate_tv_sf_data.R")
source("functions/generate_tv_sf_data.R")
source("functions/smooth_generate_tv_sf_data.R")
source("functions/loggle.graph.R")
source("functions/generate_tv_er_data.R")

tr = function(x) sum(diag(x))

# Simulate time-varying graph data using 
# hglasso_generate_tv_sf_data

set.seed(1)

# # hglasso, huge, liu_ihler, smooth, yang_peng, er

simulator = "er"

p = 100
n = 10 # Nmb of repetitions (10, 50, 100)
N = 200 # nmb of time points

# Repeat the process SM times to take into account
# the random error caused by the simulation process

SM = 10

Results = vector("list", SM)

if(simulator == "yang_peng"){
  
  example <- loggle.graph(p = p, n = n, N = N)
  
  Theta_true = example$Omega.true
  
  Sigma_true = lapply(Theta_true, as.matrix)
  
  for(i in 1:N){
    
    # This is the Cholesky decomposition of the 
    # precision matrix
    
    Sigma_true[[i]] = Cholesky(Theta_true[[i]])
    
  }
    
}

if(simulator == "smooth"){
  
  init_nodes = 4
  e_add_del = 10 # how many edges are added/removed at the change point
  cp_step = 10 # every cp_step time step, e_add_del nodes are deleted and added
  cp = seq(cp_step, N - cp_step, by = cp_step) # changepoints
  
  example = smooth_generate_tv_sf_data(n = n, 
                                       N = N, 
                                       p = p,
                                       init_nodes = init_nodes,
                                       e_add_del = e_add_del,
                                       cp_step = cp_step,
                                       cp = cp)
  
  names(example) = c("X", "A", "Omega.true", "Sigma.true")
  
}

if(simulator == "liu_ihler"){
  
  example <- generate_tv_sf_data(p = p, 
                                 n = n, 
                                 N = N, 
                                 alpha = 0.2)
}

if(simulator == "huge"){
  
  example <- huge_generate_tv_sf_data(p = p, 
                                      n = n, 
                                      N = N, 
                                      alpha = 0.2)
  
}

if(simulator == "er"){
  
  e_add_del = 10 # how many edges are added/removed at the change point
  cp_step = 10 # every cp_step time step, e_add_del nodes are deleted and added
  cp = seq(cp_step, N - cp_step, by = cp_step) # changepoints
  
  example = generate_tv_er_data(n = n, 
                                N = N, 
                                p = p,
                                e_add_del = e_add_del,
                                cp_step = cp_step,
                                cp = cp)
  
}

if(simulator == "hglasso"){
  
  example <- hglasso_generate_tv_sf_data(p = p, 
                                         n = n, 
                                         N = N, 
                                         hubnumber = 5,
                                         hubsparsity = 0.95,
                                         sparsity = 0.99,
                                         alpha = 0.01) 
}

Y = example$X

Theta_true = example$Omega.true

A_true = lapply(Theta_true, as.matrix)

for(i in 1:N){
  
  A_temp = A_true[[i]]
  
  A_temp[A_temp != 0] = 1
  
  A_true[[i]] = A_temp
  
}

for(sim_i in 1:SM){
  
  Y = vector("list", N)
  
  for(i in 1:N){
    
    if(simulator != "yang_peng"){
      Y[[i]] = t(MASS::mvrnorm(n, mu = rep(0, p), Sigma = example$Sigma.true[[i]])) 
    }
    
    if(simulator == "yang_peng"){
      Y[[i]] <- t(rmvn.sparse(n = n, mu = rep(0, p), 
                              Sigma_true[[i]]))
    }
    
  }
  
  Y = do.call(cbind, Y)
  
  Y = t(Y)
  
  S = cor(Y)
  
  #Stationary glasso, pooled data (N*n "iid" observations) 
  
  nlambda = 25
  
  R_gl_pooled = matrix(0, 12, 1)
  
  lambda.min.ratio = 0.1
  
  for(i in 1:N){
    
    lambda.max = max(max(S - diag(p)), 
                     -min(S - diag(p)))
    lambda.min = lambda.min.ratio*lambda.max
    lambdas = exp(seq(log(lambda.max), 
                      log(lambda.min), 
                      length = nlambda))
    
    res_gl = list(icov = list(), 
                  df = c())
    
    j = 1
    
    for(lambda in lambdas){
      res_gl$icov[[j]] = glassoFast::glassoFast(S = S,
                                                rho = lambda)$wi
      res_gl$df[j] = (sum(res_gl$icov[[j]] != 0) - ncol(res_gl$icov[[j]]))/2
      
      j = j + 1
      
    }
    
    loglik = c()
    
    for(j in 1:nlambda){
      loglik[j] = log(det(res_gl$icov[[j]])) - tr(S%*%res_gl$icov[[j]]) 
    }
    
    gamma = 0.5
    
    ebic = -N*loglik + log(N)*res_gl$df + 4*gamma*log(p)*res_gl$df
    
    ind = which.min(ebic)
    
    A_gl = res_gl$icov[[ind]]
    
    A_gl[A_gl != 0] = 1
    
    diag(A_gl) = 0
    
    R_gl_pooled = cbind(R_gl_pooled,
                 unlist(binary_class(A = A_true[[i]], 
                                     A_est = A_gl))
    )
    
  }
  
  R_gl_pooled = R_gl_pooled[, -1]
  
  m = rep("glasso",
          each = N)
  
  D = data.frame(TP = R_gl_pooled["TP", ],
                 FP = R_gl_pooled["FP", ],
                 FN = R_gl_pooled["FN", ],
                 TN = R_gl_pooled["TN", ],
                 Pre = R_gl_pooled["Precision", ],
                 TPR = R_gl_pooled["TPR", ],
                 MCC = R_gl_pooled["MCC", ],
                 FDR = R_gl_pooled["FDR", ],
                 FPR = R_gl_pooled["FPR", ],
                 F1 = R_gl_pooled["F1", ],
                 JI = R_gl_pooled["JI", ],
                 ED = R_gl_pooled["ED", ],
                 Method = m,
                 Sim_round = sim_i)
  
  
  Results[[sim_i]] = D
  
  cat("\r", sim_i)
  
}

Sim_res = dplyr::bind_rows(Results)

Sim_res = tidyr::pivot_longer(Sim_res, 
                              cols = 1:12, 
                              values_to = "Value",
                              names_to = "Metric")

Sim_res$N = N
Sim_res$n = n
Sim_res$p = p

##

f_name = paste0("simulations/simulations_with_more_rounds/glasso_pooled/results/", 
simulator,"/", simulator, "_mod_", n, "_", p, "_", N, ".txt")

fn = paste0(f_name)
write.table(Sim_res, file = fn, row.names = FALSE)