library(tidyverse)
library(igraph)
set.seed(12345)
p = 100 # solmujen lkm
n = 10 # aikapisteiden lkm
rho = 0.25
w1 = 0.1
w2 = 0.3
e_add_del = 5 # montako solmua poistetaan/lisataan muutospisteen jalkeen
cp = c(5) # aikapiste/pisteet, jolloin muutos tapahtuu

n_n <- diff(c(cp, n)) + 1

a <- runif(p, w1, w2)

I <- diag(1, p)

Theta <- rho*I

G <- igraph::graph_from_adjacency_matrix(matrix(1, p, p), 
                                         diag = FALSE, 
                                         mode = "undirected")

E <- igraph::as_edgelist(G)

ind <- sample(1:nrow(E), p)

Es <- E[ind, ]

Theta[Es] <- Theta[Es] - a

Theta[Es[, 2:1]] <- Theta[Es[, 2:1]] - a

b <- diag(diag(Theta), p)

diag(Theta) <- diag(Theta) + abs(rowSums(Theta - b))

G_Theta = igraph::graph_from_adjacency_matrix(Theta, diag = FALSE,
                                              mode = "undirected",
                                              weighted = TRUE)

E = igraph::as_data_frame(G_Theta)

colnames(E) <- c("from", "to", "t1")

E = E %>% 
  as_tibble() %>% 
  complete(from = 1:p, to = 1:p, fill = list(t1 = 0)) %>% 
  filter(from != to)

E_temp = matrix(0, nrow = nrow(E), ncol = n - 1)

colnames(E_temp) = paste0("t", 2:n)

E = cbind(E, E_temp)
rm(E_temp)

E[3:(cp[1] + 2)] = E[, 3]

cp = c(cp, n)

for(k in 1:c(length(cp) - 1)){
  
  add_ind = which(E[, cp[k] - 1] == 0)
  
  add_edges <- sample(add_ind, e_add_del)
  
  a_target_add <- runif(e_add_del, min = w1, max = w2)
  
  a_target_add <- matrix(0, n_n[k], e_add_del)
  
  a_target_add[1, ] <- runif(e_add_del, min = w1, max = w2)
  
  a_target_add <- apply(a_target_add, 2, 
                        function(x) seq(0, x[1], length.out = n_n[k]))
  
  a_target_add = t(a_target_add)
  
  rownames(a_target_add) = add_edges
  
  E[add_edges, c(cp[k]:cp[k + 1] + 2)] = -a_target_add
  
  E[-add_edges, c(cp[k]:cp[k + 1] + 2)] = E[-add_edges, cp[k] + 2]
  
  ##
  
  del_ind = which(E[, cp[k] - 1] != 0)
  del_edges <- sample(del_ind, e_add_del)
  
  a_target_del <- matrix(0, n_n[k], e_add_del)
  
  a_target_del[1, ] <- E[del_edges, cp[k]]
  
  a_target_del <- apply(a_target_del, 2, 
                        function(x) seq(x[1], 0, length.out = n_n[k]))
  
  a_target_del = t(a_target_del)
  
  rownames(a_target_del) = del_edges
  
  E[del_edges, c(cp[k]:cp[k + 1] + 2)] = a_target_del
  
  E[-c(del_edges, add_edges), c(cp[k]:cp[k + 1] + 2)] = E[-c(del_edges, add_edges), cp[k] + 2]
  
  
}


n_replicas = 10

ts_sim_data_100p_125n <- list()
ts_true_data_100p_125n <- list()
ts_true_sigma_100p_125n <- list()


for (j in 1:n_replicas) {
  
  ts_sim_data_100p_125n[[j]] <- list()
  ts_true_data_100p_125n[[j]] <- list()
  ts_true_sigma_100p_125n[[j]] <- list()
  
  
  for(i in 1:n){
    
    E_temp = as.matrix(E[, c(1, 2, i + 2)])
    
    G <- igraph::graph_from_edgelist(E_temp[, c(1, 2)])
    
    E(G)$weight = E_temp[, 3]
    
    Theta = igraph::as_adjacency_matrix(G, attr="weight")
    
    Theta = as.matrix(Theta)
    
    Theta = Theta + t(Theta)
    
    diag(Theta) = rho
    
    adj <- Theta
    diag(adj) = 0
    adj[adj !=0] <- 1
    ts_true_data_100p_125n[[j]][[i]] = adj
    
    b <- diag(diag(Theta), p)
    
    diag(Theta) <- diag(Theta) + abs(rowSums(Theta - b))
    
    
    Sigma = solve(Theta)
    
    ts_sim_data_100p_125n[[j]][[i]] = t(MASS::mvrnorm(n = 125, mu = rep(0, p), Sigma = cov2cor(Sigma) ))
    ts_true_sigma_100p_125n[[j]][[i]] = cov2cor(Sigma)
    
  }
}

