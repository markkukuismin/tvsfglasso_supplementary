
library(dplyr)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121160

# Whole embryos of fruit fly measured at 14 time points
# (0h, 1h, 2h, 3h, 4h, 5h, 6h, 8h, 10h, 12h, 14h, 16h, 
# 18h, 20h). Each sample was measured in biological 
# quadruplicates

# gzfiles = list.files(path = "real_data_examples/becker_fruitfly/data/GSE121160_RAW/",
#                      pattern = "*.gz", 
#                      full.names = TRUE)
# 
# rawdatapath = "real_data_examples/becker_fruitfly/data/GSE121160_RAW"
# 
# # untrar files. Note! This will deleted original files
# 
# for(i in gzfiles){
#   
#   R.utils::gunzip(i)
#   
# } 

# Load data into R

tsv_files = list.files(path = "real_data_examples/becker_fruitfly/data/GSE121160_RAW/",
                   pattern = "*.tsv", 
                   full.names = TRUE)

data = read.table(tsv_files[1], header = FALSE)

TFs = read.table("real_data_examples/becker_fruitfly/data/TFs_996.txt")

colnames(TFs) = c("ID", "name")

timepoints = rep(c("0h", "1h", "2h", "3h", "4h", "5h", 
"6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h"), each = 4)

a = c(".1", ".2", ".3", ".4")

timepoints = paste0(timepoints, a)

colnames(data) = c("gene", timepoints[1])

for(i in 2:length(tsv_files)){
  
  D = read.table(tsv_files[i], header = FALSE)

  colnames(D) = c("gene", timepoints[i])
    
  data = data %>%
    left_join(D, by = "gene")
  
}

object.size(data)/10^6

data = data[-tail(1:nrow(data), 5), ]

X = as.matrix(data[, -1])

rownames(X) = data$gene

# Choose TFs,

X = X[TFs$ID, ]
data = data[data$gene %in% TFs$ID, ]

ind = X != 0

ind = rowSums(ind)

ind = which(ind >= 20)

data = data[ind, ]

X = as.matrix(data[, -1])

ind = apply(X, 1, mean)

#ind = which(ind > 100)
ind = which(ind > 10)

data = data[ind, ]

nrow(data)

data[, -1] = log2(data[,-1] + 1)

data[, -1] = apply(data[, -1], 2, scale)

hist(unlist(data[, 2]))
hist(unlist(data[, 3]))

hist(unlist(data[, ncol(data)]))

#

plot(unlist(data[1, -1]), type = "l")

plot(unlist(data[2, -1]), type = "l")

plot(unlist(data[7, -1]), type = "l")

plot(unlist(data[11, -1]), type = "l")

plot(unlist(data[11, -1]), type = "l")

plot(unlist(data[which(data$gene %in% "FBgn0000606"), -1]),
     type = "l",
     main = "EVE")

plot(unlist(data[which(data$gene %in% "FBgn0001077"), -1]),
     type = "l",
     main = "FTZ")

plot(unlist(data[which(data$gene %in% "FBgn0001078"), -1]),
     type = "l",
     main = "FTZ-F1")

#

X = as.matrix(data[, -1])

rownames(X) = data$gene

v = apply(X, 1, function(x) var(x, na.rm = TRUE))

v = sort(v, decreasing = TRUE)[1:400]
ind = names(v)

X = X[ind, ]

dim(X)

# FBgn0000606 EVE
# FBgn0001077 FTZ
# FBgn0001078 FTZ-F1

c("FBgn0000606", "FBgn0001077", "FBgn0001078") %in% rownames(X)
