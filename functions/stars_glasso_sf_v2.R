
stars_glasso_sf_v2 <- function(X = NULL, alpha = NULL, K = 2, stars_thr = 0.1, subsamp_ratio = NULL, rep_num = 20, verbose = TRUE){
  if(is.null(X)) stop("X is missing with no default")
  if(is.null(alpha)) stop("alpha is missing with no default")
  nalpha <- length(alpha)
  n <- nrow(X)
  p <- ncol(X)
  if(is.null(subsamp_ratio)){
    if(n > 144) subsamp_ratio <- 10*sqrt(n)/n
    if(n <= 144) subsamp_ratio <- 0.8
  }
  est <- list()
  est$merge <- list()
  for(i in 1:nalpha) est$merge[[i]] <- matrix(0, p, p)
  for(i in 1:rep_num){
    if(verbose){
      mes <- paste(c("Conducting Subsampling....in progress:", 
                     floor(100*i/rep_num), "%"), collapse = "")
      cat(mes, "\r")
      flush.console()
    }
    ind <- sample(c(1:n), floor(n*subsamp_ratio), replace = FALSE)
    for(j in 1:nalpha){
      S <- cov(X[ind, ])
      A <- sf_glasso(S=S, alpha=alpha[j], K=K)
      A <- abs(A)
      A[A < 1e-04] <- 0
      A[A != 0] <- 1
      est$merge[[j]] <- est$merge[[j]] + A 
    }
    rm(ind, A)
    gc()
  }
  if(verbose){
    mes <- "Conducting Subsampling....done.                 "
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  est$variability <- rep(0, nalpha)
  for(i in 1:nalpha){
    est$merge[[i]] <- est$merge[[i]]/rep_num
    est$variability[i] <- 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(p*(p-1))
  }
  est$opt.index <- max(which.max(est$variability >= 
                                   stars_thr)[1] - 1, 1)
  return(est)
}