similarity_plots <- function(adjecency_T,adjecency_E){
  
  par(mfrow=c(2,2))
  
  samat <- adjecency_T*adjecency_E
  
  network <- graph_from_adjacency_matrix(samat, mode = "undirected",diag = F)
  
  coords<-layout.fruchterman.reingold(network)
  
  plot(network,layout=coords,main="TRUE POSITIVES", vertex.label.color="black", edge.width=3,vertex.label.dist=0.5,vertex.size=3)
  
  fp <- adjecency_E - samat 
  
  network <- graph_from_adjacency_matrix(fp, mode = "undirected",diag = F)
  
  coords<-layout.fruchterman.reingold(network)
  
  # plot it
  plot(network,layout=coords,main="FALSE POSITIVES", vertex.label.color="black", edge.width=3,vertex.label.dist=0.5,vertex.size=3)
  
  fn <- adjecency_T - samat 
  
  network <- graph_from_adjacency_matrix(fn, mode = "undirected",diag = F)
  
  coords<-layout.fruchterman.reingold(network)
  
  # plot it
  plot(network,layout=coords,main="FALSE NEGATIVES", vertex.label.color="black", edge.width=3,vertex.label.dist=0.5,vertex.size=3)
  
  
  
  network <- graph_from_adjacency_matrix(adjecency_T, mode = "undirected",diag = F)
  
  coords<-layout.fruchterman.reingold(network)
  
  # plot it
  plot(network,layout=coords,main="TRUE PLOT", vertex.label.color="black", edge.width=3,vertex.label.dist=0.5,vertex.size=3)
  
  
}

library(igraph)

graph_plot <- function(adj) {
  network <- graph_from_adjacency_matrix(adj, mode = "undirected",diag = F)
  
  coords<-layout.fruchterman.reingold(network)
  
  plot(network,layout=coords,main="", vertex.label.color="black", edge.width=3,vertex.label.dist=0.5,vertex.size=3)
  
}


similarity_stats <- function(adjecency_T,adjecency_E){
  
  samat <- adjecency_T*adjecency_E
  
  fp <- adjecency_E - samat 
  
  fn <- adjecency_T - samat 
  
  unit <- matrix(1,nrow = dim(adjecency_T)[1],ncol = dim(adjecency_T)[1])
  diag(unit) <- 0
  tn <- ( unit - adjecency_T)*(unit - adjecency_E)
  
  return( matrix(c(sum(samat)/2,sum(fn)/2,sum(fp)/2,sum(tn)/2),2,2,byrow=TRUE))
}





calculate_scores <- function(cm) {
  tp <- cm[1,1]
  tn <- cm[2,2]
  fp <- cm[2,1]
  fn <- cm[1,2]
  tpr <- tp / (tp + fn)
  tnr <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  fnr <- 1 - tpr
  fpr <- 1 - tnr
  fdr <- 1 - ppv
  FOR <- 1 - npv
  lr_plus <- tpr / fpr
  lr_neg <- fnr / tnr
  pt <- sqrt(fpr) / (sqrt(tpr) + sqrt(fpr))
  ts <- tp / (tp + fn + fp)
  acc <- (tp + tn) / (tp + tn + fn + fp)
  bal_acc <- (tpr + tnr) / 2
  F1_score <- 2 * (ppv * tpr) / (ppv + tpr)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  results <- data.frame(ACC = acc, ACC_bal = bal_acc, MCC = mcc, F1 = F1_score,
                        TPR = tpr, TNR = tnr, PPV = ppv, NPV = npv, FNR = fnr,
                        FPR = fpr, FOR = FOR, LRp = lr_plus, LRn = lr_neg)
  return(results)
}

ts_scores <- function(true,est){
  nt <- length(est)
  
  scoreboard <- matrix(NA, nt, 13)
  colnames(scoreboard) <- colnames(calculate_scores( similarity_stats(true[[1]],est[[1]])))
  for (i in 1:nt) {
    scoreboard[i,] <- as.matrix(calculate_scores( similarity_stats(true[[i]],est[[i]])))
  }
  return(scoreboard)
}




rho_res_mat_to_adj_list <- function(rho_res_mat, N, M, n_time_points,credible_level){
  
  Al <- list()
  for (i in 1:n_time_points) {
    Al[[i]] <- matrix(0,M,M)
  }
  
  Au <- list()
  for (i in 1:n_time_points) {
    Au[[i]] <- matrix(0,M,M)
  }
  
  A <- list()
  for (i in 1:n_time_points) {
    A[[i]] <- matrix(0,M,M)
  }
  
  res_ql <- matrix(apply(rho_res_mat$draws(variables = "rho_l", format = "matrix"),2, function(x) quantile(x,probs=c(1-credible_level))),nrow = n_time_points)
  #res_ql <-  res_ql - median(res_ql)
  
  res_qu <- matrix(apply(rho_res_mat$draws(variables = "rho_l", format = "matrix"),2,function(x) quantile(x,probs=c(credible_level))),nrow=n_time_points)
  #res_qu <-  res_qu + median(res_qu)
  
  for (ti in 1:n_time_points) {
    for (i in 2:M) {
      for (j in 1:(i-1)) {
        Al[[ti]][i,j] <- res_ql[ti,(j+(i-1)*(i-2)/2)]
        Al[[ti]][j,i] <- Al[[ti]][i,j]
      }  
    } 
  }
  
  for (ti in 1:n_time_points) {
    for (i in 2:M) {
      for (j in 1:(i-1)) {
        Au[[ti]][i,j] <- res_qu[ti,(j+(i-1)*(i-2)/2)]
        Au[[ti]][j,i] <- Au[[ti]][i,j]
      }  
    } 
  }
  
  for (ti in 1:n_time_points) {
    A[[ti]] <- rajat(Al[[ti]],Au[[ti]])
  }
  
  return(A)
}




rajat <- function(yla,ala){
  arvo <- yla-ala;
  i <- dim(arvo)[1]
  
  for (p in 1:i) {
    for (h in 1:i) {
      if (max(abs(yla[p,h]),abs(ala[p,h]))>abs(arvo[p,h])){
        arvo[p,h] <- 1
      }
      else{
        arvo[p,h] <- 0
      }
    }
  }
  return(arvo)
}


omega2adjancy <- function(res, rounded = 3){
  adj_list <- list()
  tt <- length(res)
  
  for (i in 1:tt){
    theta <- round(res[[i]], 3)
    theta[theta != 0] <- 1
    diag(theta) <- 0
    
    adj_list[[i]] <- theta
  }
  
  return(adj_list)
}


