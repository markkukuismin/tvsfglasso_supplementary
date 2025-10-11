
# This is the kernel estimator of the covariance matrix applied in 
# Zhou et al (2008) - Time Varying Undirected Graphs
#
# Original script copied from https://github.com/jlyang1990/loggle,
# Yang, J. & Peng, J. (2018), Estimating Time-Varying Graphical Models

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# N: number of time grid points
# pos: indices of time points used to generate sample covariance/correlation matrix
# h: bandwidth in kernel smoothed sample covariance/correlation matrix
# use_cor: if TRUE, use sample correlation matrix in model fitting,
#           if FALSE, use sample covariance matrix in model fitting

# Output ###
# S_t: list of kernel smoothed sample covariance/correlation matrices
# sdX: if use_cor = TRUE: list of standard deviations of variables
#      if use_cor = FALSE: list of 1's

kernel_cov_rep <- function(X = NULL, N = NULL, pos = NULL, h = NULL, kernel = "gaussian") {
  
  if(!is.list(X)) stop("X must be a list")
  if(is.null(X)) stop("X is missing with no default")
  if(is.null(N)) stop("N is missing with no default")
  
  if(is.null(pos)) pos <- 1:N
  
  p <- nrow(X[[1]])
  
  #if(is.null(h)) h <- 5.848/N^(1/3) # Zhou et al 2010 when changes time point 200, 200^1/3 = 5.848
  if(is.null(h)) h <- 1

  S_temp <- vector("list", N)
  
  S_t <- vector("list", length(pos))
  
  for(i in 1:N){
    
    S_temp[[i]] <- cov(t(X[[i]]))
    
  }
  
  for(i in 1:N){
    if(kernel == "epanechnikov") Kh <- pmax(3/4*(1 - ((pos - i)/((N - 1)*h))^2), 0)
    if(kernel == "gaussian") Kh <- pmax(dnorm(((pos - i)/((N - 1)*h))), 0) 
    omega <- Kh/sum(Kh)
    ind <- which(omega != 0)
    K <- matrix(0, p, p)
    for(j in ind){
      K <- K + S_temp[[pos[j]]]*omega[j]
    }
    S_t[[i]] <- K
  }
  
  result <- list(S_t = S_t)
  return(result)
}