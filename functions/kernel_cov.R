
# This is the kernel estimator of the covariance matrix applied in 
# Zhou et al (2008) - Time Varying Undirected Graphs
#
# Original script copied from https://github.com/jlyang1990/loggle,
# Yang, J. & Peng, J. (2018), Estimating Time-Varying Graphical Models

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# pos: indices of time points used to generate sample covariance/correlation matrix
# h: bandwidth in kernel smoothed sample covariance/correlation matrix
# use_cor: if TRUE, use sample correlation matrix in model fitting,
#           if FALSE, use sample covariance matrix in model fitting

# Output ###
# S_t: list of kernel smoothed sample covariance/correlation matrices
# sdX: if use_cor = TRUE: list of standard deviations of variables
#      if use_cor = FALSE: list of 1's

kernel_cov <- function(X = NULL, pos = NULL, h = NULL, use_cor = FALSE, kernel = "gaussian") {
  
  if(!is.numeric(X)) stop("X must be numeric")
  if(is.null(X)) stop("X is missing with no default")
  
  p <- nrow(X)
  N <- ncol(X)
  
  if(is.null(pos)) pos <- 1:N
  
  if(is.null(h)) h <- 5.848/N^(1/3) # (200)^(1/3) = 5.848, See Zhou et al (2010) examples
  
  sdX <- rep(1, p)
  
  if(use_cor){
    for(i in 1:p){
      sdX[i] <- sd(X[i, ])
      X[i, ] <- X[i, ]/sd(X[i, ])
    }
  }
  
  S_t <- array(0, c(p, p, N))
  
  for(i in 1:N){
    if(kernel == "epanechnikov") Kh <- pmax(3/4*(1 - ((pos - i)/((N - 1)*h))^2), 0)
    if(kernel == "gaussian") Kh <- pmax(dnorm(((pos - i)/((N - 1)*h))), 0) 
    omega <- Kh/sum(Kh)
    ind <- which(omega != 0)
    X_pos <- X[, pos[ind]]
    S_t[, , i] <- (X_pos*rep(omega[ind], rep(p, length(ind))))%*%t(X_pos)
  }
  
  result <- list(S_t = S_t, sdX = sdX)
  return(result)
}