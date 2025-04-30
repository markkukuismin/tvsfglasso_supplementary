
simulate_sf <- function(p = 100, n = 200, eta = 1.1){
  
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
  
  Theta <- Lambda2%*%L%*%Lambda2

  Sigma <- solve(Theta)
  
  mu <- rep(0, p)
  
  X <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma)
  
  return(list(G = G, Theta = Theta, Sigma = Sigma, X = X))
              
}