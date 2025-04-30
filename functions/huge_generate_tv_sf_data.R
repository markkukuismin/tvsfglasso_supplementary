# function to generate time-varying graphs and 
# observations
huge_generate_tv_sf_data <- function(p, n = 1, N, v = 0.3, u = 0.1, alpha = 0) {
  
  data <- huge::huge.generator(n = 10,
                               d = p,
                               graph = "scale-free",
                               v = v,
                               u = u,
                               verbose = FALSE) 
  
  Omega_t <- data$omega
    
  Omega.true <- vector("list", N)
  Sigma.true <- vector("list", N)
  
  A <- as.matrix(data$theta)
  
  G_t <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  
  G <- G_t
  
  ebc <- edge_betweenness(G_t, directed = FALSE)*(2/(p*(p - 1)))
  
  E <- as_edgelist(G_t)
  
  E <- E[order(E[, 1]), ]
  
  E <- cbind(E, ebc)
  
  A_t <- A
  
  B <- A_t
  
  S <- c()
  
  for(i in 1:nrow(E)){
    
    B[E[i, 1], E[i, 2]] <- B[E[i, 2], E[i, 1]] <- E[i, "ebc"]
    
  }
  
  X <- vector("list", N)
  
  for(i in 1:N){
    
    R <- matrix(runif(p*p, 
                      min = alpha, 
                      max = max(ebc)), p, p)
    
    R <- (R + t(R))/2
    
    Ind <- R < B
    
    Ind <- Ind*1
    
    A_t <- A_t - Ind
    
    A_t[A_t != 0] <- 1
    
    Omega_t <- A_t*v
    
    diag(Omega_t) <- abs(min(eigen(Omega_t)$values)) + 0.1 + u
    
    Sigma <- cov2cor(solve(Omega_t))
    
    Omega.true[[i]] <- Omega_t
    
    Sigma.true[[i]] <- Sigma
    
    X[[i]] <- t(MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
    
    S[i] = isSymmetric(Omega_t)
    
  }
  
  result <- list(Omega.true = Omega.true, X = X, S = S, Sigma.true = Sigma.true)
  return(result)
}