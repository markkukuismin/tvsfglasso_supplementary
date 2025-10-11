
K <- function(s, t, kernel = "gaussian"){
  if(kernel == "gaussian") w <- (1/sqrt(2*pi))*exp(-0.5*(abs(s - t)/h)^2)
  if(kernel == "epanechnikov") w <- (3/4)*(1 - (abs(s - t)/h)^2)
  w
}

# Y = list of length N of p x n data matrices
# s = normalized timevector (between 0 and 1)
# h = bandwidth

weighted_cov <- function(Y, s = NULL, h = NULL){
  
  p <- nrow(Y[[1]])
  
  N <- length(Y)
  
  if(is.null(h)) h <<- 5.848/N^(1/3)
  
  S_temp <- vector("list", N)
  
  S_t <- vector("list", N)
  
  for(i in 1:N){
    
    S_temp[[i]] <- cov(t(Y[[i]]))
    
  }
  
  i <- 1
  
  for(t in s){
    w = K(s, t, kernel = "epanechnikov")
    w = w/sum(w)
    Sw <- matrix(0, p, p)
    ind <- which(w != 0)
    for(j in ind){
      Sw <- Sw + S_temp[[j]]*w[j]
    }
    S_t[[i]] <- Sw
    i <- i + 1
  }
  
  S_t
  
}
