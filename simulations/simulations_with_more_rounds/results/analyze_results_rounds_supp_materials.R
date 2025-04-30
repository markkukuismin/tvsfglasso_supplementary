
library(tidyverse)
library(dplyr)

# yang_peng, huge, liu_ihler, hglasso, or smooth

simulator = "smooth"

if(simulator == "liu_ihler"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/liu_ihler", pattern = "*.txt", full.names = TRUE)
  files_temp_glasso_pooled = list.files(path = "simulations/simulations_with_more_rounds/glasso_pooled/results/liu_ihler", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "yang_peng"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/yang_peng", pattern = "*.txt", full.names = TRUE)
  files_temp_glasso_pooled = list.files(path = "simulations/simulations_with_more_rounds/glasso_pooled/results/yang_peng", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "huge"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/huge", pattern = "*.txt", full.names = TRUE)
  files_temp_glasso_pooled = list.files(path = "simulations/simulations_with_more_rounds/glasso_pooled/results/huge", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "hglasso"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/hglasso", pattern = "*.txt", full.names = TRUE)
  files_temp_glasso_pooled = list.files(path = "simulations/simulations_with_more_rounds/glasso_pooled/results/hglasso", pattern = "*.txt", full.names = TRUE)
}

if(simulator == "smooth"){
  files_temp = list.files(path = "simulations/simulations_with_more_rounds/results/smooth", pattern = "*.txt", full.names = TRUE)
  files_temp_glasso_pooled = list.files(path = "simulations/simulations_with_more_rounds/glasso_pooled/results/smooth", pattern = "*.txt", full.names = TRUE)
}

Data = lapply(files_temp, 
              function(x) read.table(x, 
                                     header = TRUE, 
                                     stringsAsFactors = TRUE))

Data = dplyr::bind_rows(Data)

#

Datagl = lapply(files_temp_glasso_pooled, 
                function(x) read.table(x, 
                                     header = TRUE, 
                                     stringsAsFactors = TRUE))

Datagl = dplyr::bind_rows(Datagl)

Datagl$Method = "pglasso"

Data = dplyr::bind_rows(Data, Datagl)

#

Data$N = as.factor(Data$N)
Data$n = as.factor(Data$n)
Data$p = as.factor(Data$p)

Data$Metric = as.factor(Data$Metric)

Data = Data %>%
  filter(Metric %in% c("FPR", "MCC", "FDR", "Pre", "TPR", "F1")) %>%
  filter(Method %in% c("pglasso", "tvglasso", "tvsfglasso"))


boxp = ggplot(data = Data, aes(x = Method, y = Value, color = n)) +
  geom_boxplot() +
  theme(legend.position = "none" ) +
  facet_wrap(~Metric, ncol = 2, scales = "free")

#boxp

Data$Metric = droplevels(Data$Metric)

boxp = ggplot(data = Data, aes(x = Method, y = Value, fill = n)) +
  geom_boxplot() +
  facet_wrap(~Metric, ncol = 2, scales = "free") +
  # Colorblind-friendly colors
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  labs(fill = expression(n[k]))

#boxp

Data_temp = Data %>%
  filter(Method %in% c("tvglasso", "tvsfglasso"))

boxp_temp = ggplot(data = Data_temp, 
                   aes(x = Method, y = Value, fill = n)) +
  geom_boxplot() +
  facet_wrap(~Metric, ncol = 2, scales = "free") +
  # Colorblind-friendly colors
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + 
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

#boxp_temp

# if(simulator == "liu_ihler"){
#   ggsave("simulations/simulations_with_more_rounds/results/figures_supplementary/liu_ihler_simulator.png",
#          plot = boxp,
#          dpi = 600,
#          width = 10,
#          height = 7)
# }

if(simulator == "liu_ihler"){
  ggsave("simulations/simulations_with_more_rounds/results/figures_supplementary/liu_ihler_simulator.pdf",
         plot = boxp,
         width = 10,
         height = 7)
}

if(simulator == "yang_peng"){
  ggsave("simulations/simulations_with_more_rounds/results/figures_supplementary/yang_peng_simulator.pdf",
         plot = boxp,
         width = 10,
         height = 7)
}

if(simulator == "huge"){
  ggsave("simulations/simulations_with_more_rounds/results/figures_supplementary/huge_simulator.pdf",
         plot = boxp,
         width = 10,
         height = 7)
}

if(simulator == "hglasso"){
  ggsave("simulations/simulations_with_more_rounds/results/figures_supplementary/hglasso_simulator.pdf",
         plot = boxp,
         width = 10,
         height = 7)
}


if(simulator == "smooth"){
  ggsave("simulations/simulations_with_more_rounds/results/figures_supplementary/smooth_simulator.pdf",
         plot = boxp,
         width = 10,
         height = 7)
}