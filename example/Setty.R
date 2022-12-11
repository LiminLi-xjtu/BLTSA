library(FNN)
library(Matrix)
library(igraph)
library(ggplot2)
library(RSpectra)

# use Dimensionality.R to reduce the dimensionality

data <- R.matlab::readMat('sim4.mat')
# data after dimensionality reduction
X <- data$DiffusionMap
X <- data.frame(dim1=X[,1],dim2=X[,2])

# calculate labels
result<-mylabelling_fix(X,k=1000)
label<-result[["label"]]
cutoff<-result[["cutoff"]]
root<-result[["root"]]

# calculate pseudotime
pseudotime = align(X,label,root=278,c(40,400,cutoff))
