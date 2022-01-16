BiocManager::install("destiny")
library(destiny)


X <- read.csv('simulation.csv',row.names= 1,header=TRUE)
D = DiffusionMap(X)
X1 <- data.frame(dim1=D@eigenvectors[,1],dim2=D@eigenvectors[,2])
write.csv(X1, "sim.csv")