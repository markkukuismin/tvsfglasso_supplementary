
glasso_sf <- function(X, alpha = NULL, K = 5){
  
  if(is.null(alpha)) alpha <- 0.1
  
  k <- 1
  
  p <- ncol(X)
  
  Theta_n <- diag(1, p)
  
  lambdaij <- matrix(0, p, p)
    
  while(k <= K){
    
    ep <- diag(Theta_n)
    
    for(i in 1:p){
      for(j in 1:p){
        
        lambdaij[i, j] <- alpha*(1/(sum(abs(Theta_n[i, -i])) + ep[i]) + 1/(sum(abs(Theta_n[j, -j])) + ep[j]))
        
      }
    }
    
    Theta <- CVXR::Variable(p, p, PSD = TRUE)
    
    if(!isSymmetric(X)) S <- cov(X)
    
    if(isSymmetric(X)) S <- X
    
    beta <- 2*alpha/ep
    
    obj <- CVXR::log_det(Theta) - CVXR::matrix_trace(S%*%Theta) - sum(lambdaij[lower.tri(lambdaij)]*abs(Theta[lower.tri(Theta)]))*2 - sum(abs(beta*diag(Theta)))
    
    prob <- CVXR::Problem(Maximize(obj))
    
    result <- CVXR::solve(prob, solver = "SCS")
    
    Theta_n <- result$getValue(Theta)
    
    k <- k + 1
    
  }
  
  return(list(res = result, Theta_hat = Theta_n))
  
}