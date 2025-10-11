
library(tidyverse)
library(dplyr)

# hglasso, huge, liu_ihler, yang_peng, smooth, er

simulator = "liu_ihler"

n = 10 # 10, 50 or 100

if(simulator == "liu_ihler"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/ebic_sensitivity", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "yang_peng"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/ebic_sensitivity", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "huge"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/ebic_sensitivity", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "hglasso"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/ebic_sensitivity", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "smooth"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/ebic_sensitivity", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "er"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/ebic_sensitivity", pattern = "*.txt", full.names = TRUE)
}

Data = lapply(files_temp, 
              function(x) read.table(x, 
                                     header = TRUE, 
                                     stringsAsFactors = TRUE))

Data = dplyr::bind_rows(Data)

Data$N = as.factor(Data$N)
Data$n = as.factor(Data$n)
Data$p = as.factor(Data$p)

ebic_gamma = c(0.1, 0.2, 0.5, 0.7, 0.9)

Data$gamma = rep(ebic_gamma, each = 200*10*12)
Data$gamma = as.factor(Data$gamma)

Data$Metric = as.factor(Data$Metric)

Data = Data %>%
  filter(Metric %in% c("FPR", "MCC", "FDR", "Pre", "TPR", "F1", "JI", "ED"))

Data$Metric = droplevels(Data$Metric)

boxp = ggplot(data = Data, aes(x = Method, y = Value, fill = gamma)) +
  geom_boxplot() +
  facet_wrap(~Metric, ncol = 2, scales = "free") +
  # Colorblind-friendly colors
  #scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  labs(fill = expression(gamma), x = NULL)

#boxp

if(simulator == "liu_ihler" & n == 10){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/liu_ihler_simulator_ebic_n10.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "liu_ihler" & n == 50){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/liu_ihler_simulator_ebic_n50.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "liu_ihler" & n == 100){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/liu_ihler_simulator_ebic_n100.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "yang_peng"){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/yang_peng_simulator_ebic.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "huge"){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/huge_simulator_ebic.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "hglasso"){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/hglasso_simulator_ebic.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}


if(simulator == "smooth"){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/smooth_simulator_ebic.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}

if(simulator == "er"){
  ggsave("simulations/simulations_with_more_rounds/ebic_sensitivity/figures/er_simulator_ebic.png",
         plot = boxp,
         dpi = 600,
         width = 10,
         height = 7)
}
