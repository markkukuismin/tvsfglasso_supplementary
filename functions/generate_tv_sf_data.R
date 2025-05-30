# function to generate time-varying graphs and 
# observations (Liu & Ihler 2011)
generate_tv_sf_data <- function(p, n = 1, N, eta = 1.1, alpha = 0) {
  
  G <- matrix(0, 4, 4)
  
  G[1, 2] <- G[2, 3] <- G[3, 4] <- 1
  
  G <- G + t(G)
  
  d <- colSums(G)
  
  prob <- d/sum(d)
  
  P <- ncol(G)
  
  while(P < p){
    
    for(i in 1:P){
      
      p_temp <- stats::rbinom(1, 1, prob = prob[i])
      
      if(p_temp == 1){
        
        ne <- rep(0, P + 1)
        
        ne[i] <- p_temp
        
        G_temp <- matrix(0, P + 1, P + 1)
        
        G_temp[1:P, 1:P] <- G
        
        G_temp[P + 1, ] <- ne
        
        G_temp[, P + 1] <- ne
        
        G <- G_temp
        
        rm(G_temp)
        
        d <- colSums(G)
        
        prob <- d/sum(d)  
        
      }
      
      P <- ncol(G)
      
    }
    
  }
  
  if(P > p) G <- G[1:p, 1:p]
  
  d <- colSums(G)
  
  D <- diag(d)
  
  L <- eta*D - G
  
  Linv <- solve(L)
  
  Lambda2 <- diag(sqrt(diag(Linv)))
  
  Omega_t <- Lambda2%*%L%*%Lambda2
  
  Omega.true <- vector("list", N)
  Sigma.true <- vector("list", N)
  
  A <- G
  
  G_t <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  
  G <- G_t
  
  ebc <- edge_betweenness(G_t, directed = FALSE)*(2/(p*(p - 1)))
  
  E <- as_edgelist(G_t)
  
  E <- E[order(E[, 1]), ]
  
  E <- cbind(E, ebc)
  
  A_t <- A
  
  B <- A_t
  
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
    
    d <- colSums(A_t)
    
    D <- diag(d)
    
    diag(D)[diag(D) == 0] <- 1 
    
    L <- eta*D - A_t
    
    Linv <- solve(L)
    
    Lambda2 <- diag(sqrt(diag(Linv)))
    
    Omega_t <- Lambda2%*%L%*%Lambda2
    
    Omega.true[[i]] <- Omega_t
    
    Sigma <- solve(Omega_t)
    
    Sigma.true[[i]] <- Sigma
    
    X[[i]] <- t(MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
    
  }
  
  result <- list(Omega.true = Omega.true, Sigma.true = Sigma.true, X = X)
  return(result)
}