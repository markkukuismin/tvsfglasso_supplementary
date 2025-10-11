
# The initial development version of the time varying scale-free glasso (tvsfgl)

tvsfgl <- function(X = NULL, lambda2 = 0.1, alpha = 0.1, K = 5){
  
  if(is.null(X)) stop("X is missing with no default")
  
  k <- 1
  
  Tt <- length(X)
  
  p <- ncol(X[[1]])
  
  Theta_n <- array(rep(c(diag(p)), Tt), c(p, p, Tt))
  
  obj_f_value <- rep(0, K)
  
  for(i in 1:Tt){
    
    assign(paste0("S", i), 
           stats::cov(X[[i]]))
           
  }
  
  k = 1
  
  while(k <= K){
    
    for(i in 1:Tt){
      
      assign(paste0("Theta", i), 
             CVXR::Variable(p, p, PSD = TRUE))
      
      assign(paste0("S", i), 
             cov(X[[i]]))
      
      ep <- diag(Theta_n[, , i])
      
      lambdaij <- matrix(0, p, p)
      
      for(ii in 1:p){
        for(jj in 1:p){
          
          lambdaij[ii, jj] <- alpha*(1/(sum(abs(Theta_n[ii, -ii, Tt])) + ep[ii]) + 1/(sum(abs(Theta_n[jj, -jj, Tt])) + ep[jj]))
          
        }
      }
      
      beta <- 2*alpha/ep
      
      diag(lambdaij) <- beta
      
      assign(paste0("lambdaM", i),
             lambdaij)
      
    }
    
    obj_text <- paste0("CVXR::log_det(Theta1) - CVXR::matrix_trace(S1 %*% Theta1) - sum(abs(lambdaM1[lower.tri(lambdaM1)]*Theta1[lower.tri(Theta1)]))*2")
    
    for(i in 2:Tt){
      
      obj_text <- paste0(obj_text, 
                         " + CVXR::log_det(Theta", i, ") - CVXR::matrix_trace(S", i, "%*% Theta", i, ") - sum(abs(lambdaM", i, "[lower.tri(lambdaM", i, ")]*Theta", i, "[lower.tri(Theta", i, ")]))*2 - lambda2*sum(abs(Theta", i, "[lower.tri(Theta", i, ")] - Theta", i-1, "[lower.tri(Theta", i-1, ")]))*2")
      
    }
    
    obj <- eval(parse(text = obj_text))
    
    prob <- CVXR::Problem(Maximize(obj))
    
    result <- CVXR::solve(prob, solver = "SCS", warm_start = TRUE)
    
    Theta_n <- array(0, c(p, p, Tt))
    
    for(i in 1:Tt){
      
      res <- paste0("result$getValue(Theta", i, ")")
      
      Theta_n[,,i] <- eval(parse(text = res))
      
    }
    
    obj_f_value[k] <- result$value
    
    k <- k + 1
    
  }
  
  return(list(Theta = Theta_n, value = obj_f_value))
  
}