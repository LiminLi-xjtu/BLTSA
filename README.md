# BLTSA
BLTSA (branched local tangent space alignment) is a method to infer single cell pseudotime for multi-furcation trajectories. By assuming that single cells are lying on a low-dimensional self-intersecting manifold, BLTSA first identifies the tip and branching cells in the trajectory  based on cells' local Euclidean neighborhoods, and clusters all the cells to different branches iteratively.  Local  coordinates within the tangent spaces determined by each cell's  are then computed.  The global coordinates for all the single cells  are  finally obtained by aligning the local coordinates  based on the tangent spaces. 

#### Dependencies
BLTSA is performed in R, and the following dependent packages need to be installed before using BLTSA.

     install.packages("FNN")
     install.packages("igraph")
     install.packages("RSpectra")
     BiocManager::install("destiny")
     
     library(FNN)
     library(Matrix)
     library(igraph)
     library(RSpectra)
     library(destiny)
    
#### Usage

1.The first dimension reduction part can be performed in the Dimensionalith.R file.

2.All the dependent functions including spectral clustering.R, DISTMAT.R, mylabelling_fix.R, align.R need to be imported first.

3.The parameters label,root, and cutoff can be obtained with the following code:

     result = mylabelling_fix(X)
     label = result[["label"]]
     cutoff = result[["cutoff"]]
     root = result[["root"]]
     
4.The pseudotime can be calculated as follows:

     pseudotime = align(X,label,root,c(40,100,cutoff))
