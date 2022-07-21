library(FNN)
library(Matrix)
library(igraph)
library(igraph)
library(ggplot2)

# use Dimensionality.R to reduce the dimensionality

data <- R.matlab::readMat('sim3.mat')
X <- data$sim
X <- data.frame(dim1=X[,1],dim2=X[,2])

# calculate labels
result<-mylabelling_fix(X)
label<-result[["label"]]
cutoff<-result[["cutoff"]]
root<-result[["root"]]
plot(X,col=label)

# calculate pseudotime
pseudotime = align(X,label,root,c(40,100,cutoff))
ggplot(X, aes(x = dim1, y = dim2, color = pseudotime)) + geom_point() + scale_color_viridis_c(option = "viridis")
