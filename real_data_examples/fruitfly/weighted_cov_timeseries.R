
# X = p x N data matrix
# s = normalized timevector (between 0 and 1)
# pos = time points used to generate sample covariance matrix
# h = bandwidth
# use_cor = compute sample correlation
# kernel = kernel function

weighted_cov_timeseries <- function(X = NULL, 
                                    s = NULL,
                                    pos = NULL,
                                    h = NULL,
                                    use_cor = FALSE, 
                                    kernel = "gaussian") {
  
  if(!is.numeric(X)) stop("X must be numeric")
  if(is.null(X)) stop("X is missing with no default")
  
  p <- nrow(X)
  N <- ncol(X)
  
  if(is.null(s)) stop("s is missing with no default")
  if(is.null(pos)) stop("pos is missing with no default")
  
  if(is.null(h)) h <- 5.848/N^(1/3) # (200)^(1/3) = 5.848, See Zhou et al (2010) examples
  
  sdX <- rep(1, p)
  
  if(use_cor){
    for(i in 1:p){
      sdX[i] <- sd(X[i, ])
      X[i, ] <- X[i, ]/sd(X[i, ])
    }
  }
  
  S_t <- array(0, c(p, p, length(pos)))
  
  i <- 1
  
  for(t in pos){
    if(kernel == "epanechnikov") w <- pmax(3/4*(1 - ((s - t)/h)^2), 0)
    if(kernel == "gaussian") w <- pmax(dnorm((abs(s - t)/h)), 0) 
    w <- w/sum(w)
    ind <- which(w != 0)
    X_pos <- X[, ind]
    S_t[, , i] <- (X_pos*rep(w[ind], rep(p, length(ind))))%*%t(X_pos)
    i <- i + 1
  }
  
  result <- list(S_t = S_t, sdX = sdX)
  return(result)
}
