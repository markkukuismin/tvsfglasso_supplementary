
sf_glasso <- function(S = NULL, alpha = 0.1, K = 5){
  
  if(is.null(S)) stop("S is missing with no default")
  
  k <- 1
  
  p <- ncol(S)
  
  Theta_n <- diag(1, p)
  
  while(k <= K){
    
    ep <- diag(Theta_n)
    
    row_n = rowSums(abs(Theta_n) - diag(diag(abs(Theta_n))))
    
    xr <- 1/(row_n + ep)
    lambdaM <- outer(xr, xr, FUN = "+")
    lambdaM <- alpha*lambdaM
    
    beta <- 2*alpha/ep
    
    diag(lambdaM) <- beta
    
    results <- glassoFast::glassoFast(S, rho = lambdaM)
    
    Theta_n <- results$wi
    
    k <- k + 1
    
  }
  
  Theta_n
  
}