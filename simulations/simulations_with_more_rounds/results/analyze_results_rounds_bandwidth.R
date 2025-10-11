
library(tidyverse)
library(dplyr)

# yang_peng, huge, liu_ihler, hglasso, smooth, er

simulator = "er"

n = 10 # 10, 50, 100
nr = paste0("d_", n, "_")

rmatch = c(nr, "*.txt")

if(simulator == "liu_ihler"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/liu_ihler/tvglasso_and_tvsfglasso_bandwidth_res", pattern = rmatch, full.names = TRUE)
}

if(simulator == "yang_peng"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/yang_peng/tvglasso_and_tvsfglasso_bandwidth_res", pattern = rmatch, full.names = TRUE)
}

if(simulator == "huge"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/huge/tvglasso_and_tvsfglasso_bandwidth_res", pattern = rmatch, full.names = TRUE)
}

if(simulator == "hglasso"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/hglasso/tvglasso_and_tvsfglasso_bandwidth_res", pattern = rmatch, full.names = TRUE)
}

if(simulator == "smooth"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/smooth/tvglasso_and_tvsfglasso_bandwidth_res", pattern = rmatch, full.names = TRUE)
}

if(simulator == "er"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/er/tvglasso_and_tvsfglasso_bandwidth_res", pattern = rmatch, full.names = TRUE)
}

Data = lapply(files_temp, 
              function(x) read.table(x, 
                                     header = TRUE, 
                                     stringsAsFactors = TRUE))

Data = dplyr::bind_rows(Data)

#

Data$N = as.factor(Data$N)
Data$n = as.factor(Data$n)
Data$p = as.factor(Data$p)
Data$h = as.factor(Data$h)

Data$Metric = as.factor(Data$Metric)

boxp = ggplot(data = Data, aes(x = Method, y = Value, color = h)) +
  geom_boxplot() +
  theme(legend.position = "none" ) +
  facet_wrap(~Metric, ncol = 2, scales = "free")

#boxp

Data = Data %>%
  filter(Metric %in% c("FPR", "MCC", "FDR", "Pre", "TPR", "F1", "JI", "ED"))

Data$Metric = droplevels(Data$Metric)

boxp = ggplot(data = Data, aes(x = Method, y = Value, fill = h)) +
  geom_boxplot() +
  facet_wrap(~Metric, ncol = 2, scales = "free") +
  # Colorblind-friendly colors
  #scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  labs(fill = expression(h), x = NULL)

#boxp

if(simulator == "liu_ihler"){
  fn = paste0("simulations/simulations_with_more_rounds/results/figures/bandwidth_figures/liu_ihler_simulator_", n, ".png")
  ggsave(fn,
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "yang_peng"){
  fn = paste0("simulations/simulations_with_more_rounds/results/figures/bandwidth_figures/yang_peng_simulator_", n, ".png")
  ggsave(fn,
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "huge"){
  fn = paste0("simulations/simulations_with_more_rounds/results/figures/bandwidth_figures/huge_simulator_", n,".png")
  ggsave(fn,
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "hglasso"){
  fn = paste0("simulations/simulations_with_more_rounds/results/figures/bandwidth_figures/hglasso_simulator_", n, ".png")
  ggsave(fn,
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}


if(simulator == "smooth"){
  fn = paste0("simulations/simulations_with_more_rounds/results/figures/bandwidth_figures/smooth_simulator_", n, ".png")
  ggsave(fn,
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}


if(simulator == "er"){
  fn = paste0("simulations/simulations_with_more_rounds/results/figures/bandwidth_figures/er_simulator_", n, ".png")
  ggsave(fn,
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}
