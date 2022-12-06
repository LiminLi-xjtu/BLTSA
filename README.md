# BLTSA
BLTSA (branched local tangent space alignment) is a method to infer single cell pseudotime for multi-furcation trajectories. By assuming that single cells are lying on a low-dimensional self-intersecting manifold, BLTSA first identifies the tip and branching cells in the trajectory  based on cells' local Euclidean neighborhoods, and clusters all the cells to different branches iteratively.  Local  coordinates within the tangent spaces determined by each cell's  are then computed.  The global coordinates for all the single cells  are  finally obtained by aligning the local coordinates  based on the tangent spaces. 

#### Dependencies
BLTSA is performed in R, and the following dependencies need to be installed before using BLTSA.

     install.packages("FNN")
     install.packages("igraph")
     install.packages("ggplot2")
     install.packages("RSpectra")
     BiocManager::install("destiny")
     
     library(FNN)
     library(Matrix)
     library(igraph)
     library(ggplot2)
     library(RSpectra)
     library(destiny)
    
#### Usage
1.The first dimension reduction part can be performed in  file.
