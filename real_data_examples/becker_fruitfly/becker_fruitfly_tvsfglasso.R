
library(ggplot2)
library(tidyverse)
library(igraph)
library(huge)
library(Matrix)
library(sparseMVN)
library(matrixcalc)

source("functions/kernel_cov_rep_v2.R")
source("functions/sf_glasso.R")

tr = function(x) sum(diag(x))

stages = timepoints

genes = rownames(X)

timevector = rep(c("0h", "1h", "2h", "3h", "4h", "5h", 
                   "6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h"), each = 4)

stages_uniq = unique(timevector)

N = length(stages_uniq)

X = t(X)

dim(X)

p = ncol(X)

Y = vector("list", N)

j = 1

for(s in stages_uniq){
  
  ind = stringr::str_detect(stages, paste0("^", s))
  
  Y[[j]] = t(X[ind, ])
  
  j = j + 1
  
}

lapply(Y, dim)

# time varying scale-free glasso

pos = c(0:6, seq(8, 20, by = 2))

pos
length(pos)
N

S_t = kernel_cov_rep_v2(X = Y, 
                        N = N,
                        pos = pos,
                        kernel = "epanechnikov")$S_t

nalpha = 50
alphas = seq(log(0.05), log(0.15), length.out = nalpha)
alphas = exp(alphas)

Theta_tvsfgl = array(0, c(p, p, N))

A_tvsfgl = vector("list", N)

for(i in 1:N){
  
  # Select the optimal network using eBIC,
  
  Theta_temp = array(0, c(p, p, nalpha))
  
  gamma = 0.5
  
  loglik = ebic = rep(0, nalpha)
  
  j = 1
  
  for(alpha in alphas){
    
    Theta_temp[,, j] = sf_glasso(S = S_t[[i]], 
                                 alpha = alpha,
                                 K = 2) 
    
    loglik[j] = log(det(Theta_temp[,, j])) - 
      tr(S_t[[i]]%*%Theta_temp[,, j])
    
    df = (sum(Theta_temp[,, j] != 0) - 
            ncol(Theta_temp[,, j]))/2
    
    ebic[j] = -N*loglik[j] + log(N)*df + 4*gamma*log(p)*df
    
    j = j + 1
    
  }
  
  plot(alphas, ebic)
  
  ind = which.min(ebic)
  
  Theta_tvsfgl[,,i] = Theta_temp[,, ind]
  
  A_temp = Theta_tvsfgl[,, i]
  
  A_temp = abs(A_temp)
  
  A_temp[A_temp != 0] = 1
  
  diag(A_temp) = 0
  
  A_tvsfgl[[i]] = A_temp
  
  cat("\r", i)
  
}

d_hat = d_hub = c()

Degree_all = matrix(0, N, p)

for(i in 1:N){
  
  G = graph_from_adjacency_matrix(A_tvsfgl[[i]],
                                  mode = "undirected",
                                  diag = FALSE)
  
  if(i == 1) init_h = which.max(degree(G))
  
  Degree_all[i, ] = degree(G)
  
  d_hat = c(d_hat, sum(degree(G)))
  
  d_hub = c(d_hub, degree(G)[init_h])
  
}

plot(pos, d_hat, type = "b")
abline(v = pos,
       lty = 2)

mtext(side = 3, 
      line = 0.5,
      at = pos,
      text = stages_uniq)

#

plot(pos, d_hub, type = "b")
abline(v = pos,
       lty = 2)

mtext(side = 3, 
      line = 0.5,
      at = pos,
      text = stages_uniq)


##

# fix layout,

A_temp = A_tvsfgl[[1]]

ind = which(colSums(A_temp) == 0)

A_temp = A_temp[-ind, -ind]

G = graph_from_adjacency_matrix(A_temp,
                                mode = "undirected",
                                diag = FALSE)

l = igraph::layout_with_fr(G)

plot(G, 
     vertex.label = NA, 
     vertex.size = 5,
     main = stages_uniq[1],
     layout = l)

for(i in 1:(N-1)){
  
  G_name = paste0("real_data_examples/becker_fruitfly/figures/diff_graphs/TF_net_and_diff", i, ".png")
  
  png(G_name, 
       width = 10, 
       height = 7, 
       units = 'in', 
       res = 300)
  
  A_prev = A_tvsfgl[[i]]
  #ind = which(colSums(A_prev) == 0)
  #A_prev = A_prev[-ind, -ind]
  
  A_next = A_tvsfgl[[i+1]]
  #ind = which(colSums(A_next) == 0)
  #A_next = A_next[-ind, -ind]
      
  G_prev = graph_from_adjacency_matrix(A_prev,
                                       mode = "undirected",
                                       diag = FALSE)
  
  G_next = graph_from_adjacency_matrix(A_next,
                                       mode = "undirected",
                                       diag = FALSE)
      
  diff_G = igraph::difference(G_next, G_prev)
  
  time_change = paste0(stages_uniq[i+1], "-", stages_uniq[i])
  
  par(mfrow = c(1, 3))
  
  plot(G_prev, 
        vertex.label = NA, 
        vertex.size = 5,
        main = stages_uniq[i],
        layout = l)
  
  plot(G_next, 
       vertex.label = NA, 
       vertex.size = 5,
       main = stages_uniq[i+1],
       layout = l)
  
  plot(diff_G, 
       vertex.size = 5, 
       vertex.label = NA, 
       main = time_change,
       layout = l)
    
  dev.off()
  
}

par(mfrow = c(1, 1))

##

# neighborhood of initial node with highest degree

for(i in 1:N){
  
  A_temp = A_tvsfgl[[i]]
  
  ne_hub = which(A_temp[, init_h] != 0)
  
  ne_hub = c(ne_hub, init_h)
  
  ne_hub = sort(ne_hub)
  
  A_temp = A_temp[ne_hub, ne_hub]
  
  sG = graph_from_adjacency_matrix(A_temp,
                                   mode = "undirected",
                                   diag = FALSE)
  
  plot(sG, 
       vertex.label = NA, 
       vertex.size = 10,
       main = stages_uniq[i])
  
}

##

Data = data.frame(Degree = d_hat,
                  Time = stages_uniq)

Data$Time = factor(Data$Time, levels = stages_uniq)

dplot = ggplot(data = Data, aes(Time, Degree, group = 1)) + 
  geom_line(linewidth = 1.5) +
  theme(plot.title = element_text(face = "bold")) + 
  labs(title = "(A)")

dplot

colnames(Degree_all) = colnames(X)

vd = apply(Degree_all, 2, var)

plot(vd, type = "h")

topv = order(vd, decreasing = TRUE)[1:5]

topv = names(vd)[topv]

Data = Degree_all[, topv]

Data = as.data.frame(Data)

Data$Time = stages_uniq

Data$Time = factor(Data$Time, levels = stages_uniq)

Data = pivot_longer(Data, 
                    cols = 1:5, 
                    names_to = "ID",
                    values_to = "Degree")

Data$ID = as.factor(Data$ID)

topvgplot = ggplot(data = Data, 
                   aes(x = Time, y = Degree, group = ID)) +
  geom_line(aes(linetype = ID)) + 
  theme(legend.position="bottom",
        plot.title = element_text(face = "bold")) +
  labs(title = "(B)")

topvgplot

gridplot = gridExtra::grid.arrange(dplot, topvgplot, nrow = 2)

ggsave("real_data_examples/becker_fruitfly/figures/degree_vs_time.png",
       plot = gridplot,
       dpi = 600,
       width = 10,
       height = 7)

##

lmtest = Degree_all

lmtest = as.data.frame(lmtest)

lmtest$N = 1:N

b = c()

for(i in 1:p){
  
  fo = as.formula(paste(genes[i], "~ N"))
  
  mod = lm(fo, data = lmtest)
  
  b[i] = mod$coefficients[2]
  
}

toplm = genes[order(b, decreasing = TRUE)[1:5]]

Data = Degree_all[, toplm]

Data = as.data.frame(Data)

Data$Time = stages_uniq

Data$Time = factor(Data$Time, levels = stages_uniq)

Data = pivot_longer(Data, 
                    cols = 1:5, 
                    names_to = "ID",
                    values_to = "Degree")

Data$ID = as.factor(Data$ID)

toplmplot = ggplot(data = Data, 
                   aes(x = Time, y = Degree, group = ID)) +
  geom_line(aes(linetype = ID)) + 
  theme(legend.position="bottom",
        plot.title = element_text(face = "bold")) +
  labs(title = "(B)")

toplmplot

gridplot = gridExtra::grid.arrange(dplot, toplmplot, nrow = 2)

ggsave("real_data_examples/becker_fruitfly/figures/degree_vs_time_lm.png",
       plot = gridplot,
       dpi = 600,
       width = 10,
       height = 7)

# Degree distribution of the first and the last graph,

Data = data.frame(Time = rep(stages_uniq, each = p),
                  Degree = c(t(Degree_all)))

Data$Time = as.factor(Data$Time)

Data_0h = subset(Data, Time == "0h")

ttheme = theme(plot.title = element_text(size=17),
               axis.title.y = element_text(size = 15),
               axis.title.x = element_text(size = 15),
               axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size = 12))

den_p_0h = ggplot(data = Data_0h, aes(x = Degree)) + 
  geom_density(fill = "lightblue") +
  ylab("Density") +
  labs(title = "(A)") +
  ylim(0, 0.13) +
  xlim(0, 85) +
  annotate("text", 
           x=75, y=0.12, 
           label= "0h", size = 7) +
  ttheme

den_p_0h

Data_20h = subset(Data, Time == "20h")

den_p_20h = ggplot(data = Data_20h, aes(x = Degree)) + 
  geom_density(fill = "lightblue") +
  ylab("Density") +
  labs(title = "(B)") +
  ylim(0, 0.13) +
  xlim(0, 85) +
  annotate("text", 
           x=75, y=0.12, 
           label= "20h", size = 7) +
  ttheme

den_p_20h

gridplot = gridExtra::grid.arrange(den_p_0h, 
                                   den_p_20h, 
                                   ncol = 2)

ggsave("real_data_examples/becker_fruitfly/figures/degree_density.png",
       plot = gridplot,
       dpi = 600,
       width = 10,
       height = 7)

# Difference between consecutive graphs

G_list = vector("list", N)

for(i in 1:N){
  
  G = graph_from_adjacency_matrix(A_tvsfgl[[i]],
                                  mode = "undirected",
                                  diag = FALSE)
  
  G_list[[i]] = G
  
}

ind = matrix(c(12, 11,
               13, 12,
               14, 13), ncol = 2, byrow = TRUE)

stages_uniq[12:14]

D_diff = vector("list", 3)

for(i in 1:3){
  
  g = igraph::difference(G_list[[ind[i, 1]]],
                         G_list[[ind[i, 2]]])
  
  V(g)$label = genes
  
  V(g)$color = "white"
  
  g = delete_vertices(g, degree(g) == 0)
  
  D_diff[[i]] = g
  
}


png("real_data_examples/becker_fruitfly/figures/becker_fruitfly_G_diff.png", 
     width = 10, 
     height = 7, 
     units = 'in', 
     res = 600)

par(mfrow = c(1, 3))

coords = layout_in_circle(D_diff[[1]])

plot(D_diff[[1]], 
     layout = coords,
     vertex.size = 4,
     edge.width = 4,
     vertex.label.cex = 1.5,
     vertex.label.color = "black")
title("(A)", cex.main = 3)

coords = layout_in_circle(D_diff[[2]])

plot(D_diff[[2]], 
     layout = coords,
     vertex.size = 4,
     edge.width = 4,
     vertex.label.cex = 1.5,
     vertex.label.color = "black")
title("(B)", cex.main = 3)

coords = layout_in_circle(D_diff[[3]])

plot(D_diff[[3]], 
     layout = coords,
     vertex.size = 4,
     edge.width = 4,
     vertex.label.cex = 1.5,
     vertex.label.color = "black")
title("(C)", cex.main = 3)

dev.off()
