# function to generate time-varying graphs and 
# observations
generate_tv_er_data <- function(p, n = 1, N = 10, v = 0.3, u = 0.1, rho = 0.25, w1 = 0.1, w2 = 0.3, cp_step = NULL, cp = NULL, e_add_dell = 5) {
  
  data <- huge::huge.generator(n = 10,
                               d = p,
                               graph = "random",
                               v = v,
                               u = u,
                               verbose = FALSE) 
  
  n_n <- diff(c(cp, N)) + 1
  
  Omega_t <- data$omega
    
  Omega.true <- vector("list", N)
  Sigma.true <- vector("list", N)
  
  I <- diag(1, p)
  
  Theta_temp <- rho*I
  
  A <- as.matrix(data$theta)
  
  G_t <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  
  nedges <- igraph::gsize(G_t)
  
  a <- runif(nedges, w1, w2)
  
  G <- G_t
  
  E <- igraph::as_edgelist(G_t)
  
  E <- E[order(E[,1]), ]
  
  rownames(E) <- 1:nrow(E)
  
  G <- igraph::graph_from_edgelist(E,
                                   directed = FALSE)
  
  ebc <- edge_betweenness(G, directed = FALSE)
  
  ebc <- ebc/sum(ebc)
  
  ebc <- 1 - ebc
  
  ebc <- ebc/sum(ebc)
  
  Es <- E
  
  Theta_temp[Es] <- Theta_temp[Es] - a
  
  Theta_temp[Es[, 2:1]] <- Theta_temp[Es[, 2:1]] - a
  
  b <- diag(diag(Theta_temp), p)
  
  diag(Theta_temp) <- diag(Theta_temp) + abs(rowSums(Theta_temp - b))
  
  G_Theta_temp <- igraph::graph_from_adjacency_matrix(Theta_temp, 
                                                      diag = FALSE,
                                                      mode = "undirected",
                                                      weighted = TRUE)
  
  E <- igraph::as_data_frame(G_Theta_temp)
  
  colnames(E) <- c("from", "to", "t1")
  
  E_temp <- matrix(0, nrow = nrow(E), ncol = N - 1)
  
  colnames(E_temp) <- paste0("t", 2:N)
  
  E <- cbind(E, E_temp)
  rm(E_temp)
  
  E[3:(cp[1] + 2)] <- E[, 3]
  
  cp <- c(cp, N)
  
  CM <- matrix((cp[1] + 1):N, ncol = cp_step, byrow = T)
  
  add_edges <- NULL
  
  for(k in 1:nrow(CM)){
    
    if(k >= 2){
      
      add_ind <- which(E[, CM[k, 1] + 1] == 0)
      
      add_edges <- add_ind
      
      a_target_add <- matrix(0, n_n[k], e_add_del)
      
      a_target_add[1, ] <- runif(e_add_del, 
                                 min = w1, 
                                 max = w2)
      
      a_target_add <- apply(a_target_add, 2, 
                            function(x) seq(0, x[1], length.out = n_n[k]))
      
      a_target_add <- t(a_target_add)
      
      rownames(a_target_add) <- add_edges
      
      E[add_edges, (CM[k, 1] - 1):CM[k, cp_step] + 2] <- -a_target_add
      
      E[-add_edges, (CM[k, 1] - 1):CM[k, cp_step] + 2] <- E[-add_edges, CM[k, 1] + 1]
      
      
    }
    
    ##
    
    del_ind <- which(E[, CM[k, 1] + 1] != 0)
    del_edges <- sample(del_ind, 
                        e_add_del,
                        prob = ebc[del_ind])
    
    a_target_del <- matrix(0, n_n[k], e_add_del)
    
    a_target_del[1, ] <- E[del_edges, CM[k, 1] + 1]
    
    a_target_del <- apply(a_target_del, 2, 
                          function(x) seq(x[1], 0, length.out = n_n[k]))
    
    a_target_del <- t(a_target_del)
    
    rownames(a_target_del) <- del_edges
    
    E[del_edges, (CM[k, 1] - 1):CM[k, cp_step] + 2] <- a_target_del
    
    if(!is.null(add_edges)) E[-c(del_edges, add_edges), (CM[k, 1] - 1):CM[k, cp_step] + 2] <- E[-c(del_edges, add_edges), CM[k, 1] + 1]
    if(is.null(add_edges)) E[-del_edges, (CM[k, 1] - 1):CM[k, cp_step] + 2] <- E[-del_edges, CM[k, 1] + 1]
    
  }
  
  ts_sim_data <- vector("list", N)
  ts_true_data <- vector("list", N)
  Sigma <- vector("list", N)
  Theta <- vector("list", N)
  
  for(i in 1:N){
    
    E_temp <- as.matrix(E[, c(1, 2, i + 2)])
    
    G <- igraph::graph_from_edgelist(E_temp[, c(1, 2)])
    
    E(G)$weight <- E_temp[, 3]
    
    Theta_temp <- igraph::as_adjacency_matrix(G, attr="weight")
    
    Theta_temp <- as.matrix(Theta_temp)
    
    Theta_temp <- Theta_temp + t(Theta_temp)
    
    diag(Theta_temp) <- rho
    
    adj <- Theta_temp
    diag(adj) <- 0
    adj[adj != 0] <- 1
    ts_true_data[[i]] <- adj
    
    b <- diag(diag(Theta_temp), p)
    
    diag(Theta_temp) <- diag(Theta_temp) + abs(rowSums(Theta_temp - b))
    
    Theta[[i]] <- Theta_temp
    
    Sigma_temp <- solve(Theta_temp)
    
    Sigma[[i]] <- Sigma_temp
    
    ts_sim_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cov2cor(Sigma_temp))
    
  }
  
  return(list(X = ts_sim_data, A = ts_true_data, Omega.true = Theta, Sigma.true = Sigma))
  
}
