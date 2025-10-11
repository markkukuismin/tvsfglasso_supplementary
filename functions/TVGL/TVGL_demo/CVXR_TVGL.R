
library(CVXR)

TVGL <- function(X, lambda1 = 0.1, lambda2 = 0.1, transpose = TRUE){
  
  Tt <- length(X)
  
  for(i in 1:Tt){
    
    assign(paste0("omega", i), 
           CVXR::Variable(p, p, PSD = TRUE))
    if(transpose==TRUE){
      assign(paste0("S", i), 
             cov(t(X[[i]])))
    }else{
    assign(paste0("S", i), 
           cov(X[[i]]))
    }
    
  }
  
  obj_text <- paste0("CVXR::log_det(omega1) - CVXR::matrix_trace(S1 %*% omega1) - lambda1 * sum(abs(omega1[lower.tri(omega1)])) * 2")
  
  for(i in 2:Tt){
    
    obj_text <- paste0(obj_text, 
                      " + CVXR::log_det(omega", i, ") - CVXR::matrix_trace(S", i, "%*% omega", i, ") - lambda1 * sum(abs(omega", i, "[lower.tri(omega", i, ")])) * 2 - lambda2 * sum(abs(omega", i, "[lower.tri(omega", i, ")]-omega", i-1, "[lower.tri(omega", i-1, ")])) * 2")
    
  }
  
  obj <- eval(parse(text = obj_text))
  
  prob <- CVXR::Problem(Maximize(obj))
  
  result <- CVXR::solve(prob, solver = "SCS")
  
  Omega_est <- list()
  
  for(i in 1:Tt){
    
    res <- paste0("result$getValue(omega", i, ")")
    
    Omega_est[[i]] <- eval(parse(text = res))
    
  }
  
  return(list(res = result, omega = Omega_est))
  
}

p = 10
n = 20

X = list()

for(i in 1:6){
  
  X[[i]] <- matrix(rnorm(n*p), nrow = n, ncol  = p)
  
}

test = TVGL(X, transpose = FALSE)

CVX = test$res

CVX$status

omega6 = test$omega[[6]]

norm(omega6, type = "F")

isSymmetric(omega6)

##

test2 = TVGL(X, 
             lambda1 = 0.5, 
             lambda2 = 0.5, 
             transpose = FALSE)

CVX2 = test2$res

CVX2$status

omega6.2 = test2$omega[[6]]

norm(omega6.2, type = "F")

isSymmetric(omega6.2)

norm(omega6 - omega6.2, type = "F")


######################################

test = TVGL(ts_sim_data_100p_125n[[1]],lambda1 = 0.2)
t_omega <- test$omega
rm(test)
ts_scores(ts_true_data_100p_125n[[1]], omega2adjancy(t_omega))

test2 = TVGL(ts_sim_data_100p_125n[[2]],lambda1 = 0.2)
t_omega2 <- test2$omega
rm(test2)
ts_scores(ts_true_data_100p_125n[[2]], omega2adjancy(t_omega2))

