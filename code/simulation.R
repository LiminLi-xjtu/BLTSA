library(ggplot2)
library(dplyr)
library(magrittr)

cs_sigmoid <- function(t, phi, k, delta) {
  return( 2 * phi / (1 + exp(-k*(t - delta)))) 
}

a<-function (C = 100, G = 40, p_transient = 0, zero_negative = TRUE, 
             model_dropout = FALSE, lambda = 1) 
{
  branch <- rbinom(C, 1, 0.5)
  gsd <- 2*1/rgamma(G, 2, 2)
  #gsd <- 1.5*rpois(G,3)
  k <- replicate(2, runif(G, 5, 10) * sample(c(-1, 1), G, replace = TRUE))
  phi <- replicate(2, runif(G, 5, 10))
  delta <- replicate(2, runif(G, 0.5, 1))
  inds <- 1:(G/2)
  inds2 <- (G/2 + 1):G
  k[inds, 2] <- k[inds, 1]
  k[inds2, ] <- t(apply(k[inds2, ], 1, function(r) r * sample(c(0, 
                                                                1))))
  phi[, 1] <- phi[, 2]
  delta[inds, 2] <- delta[inds, 1] <- runif(G/2, 0, 0.5)
  for (r in inds2) {
    whichzero <- which(k[r, ] == 0)
    nonzero <- which(k[r, ] != 0)
    k_sign <- sign(k[r, nonzero])
    if (k_sign == 1) {
      phi[r, whichzero] <- 0
    }
    else {
      phi[r, whichzero] <- 2 * phi[r, nonzero]
    }
  }
  pst <- runif(C)
  X <- sapply(seq_along(branch), function(i) {
    k_i <- k[, branch[i] + 1]
    phi_i <- phi[, branch[i] + 1]
    delta_i <- delta[, branch[i] + 1]
    mu <- cs_sigmoid(pst[i], phi_i, k_i, delta_i)
    rnorm(length(mu), mu, gsd)
  })
  transient_genes <- sample(C, round(p_transient * C))
  transient_genes_common <- intersect(transient_genes, inds)
  transient_genes_bifurcating <- intersect(transient_genes, 
                                           inds2)
  if (length(transient_genes_common) > 0) {
    X[transient_genes_common, ] <- t(sapply(transient_genes_common, 
                                            function(g) {
                                              scale <- rlnorm(1, log(0.05), 0.5)
                                              reverse <- sample(c(TRUE, FALSE), 1)
                                              mu <- 2 * phi[g, 1] * transient(pst, scale = scale, 
                                                                              reverse = reverse)
                                              rnorm(length(mu), mu, gsd[g])
                                            }))
  }
  if (length(transient_genes_bifurcating) > 0) {
    X[transient_genes_bifurcating, ] <- t(sapply(transient_genes_bifurcating, 
                                                 function(g) {
                                                   which_nonzero <- which(k[g, ] != 0)
                                                   scale <- rlnorm(1, log(0.05), 0.3)
                                                   reverse <- k[g, which_nonzero] < 0
                                                   mu <- 2 * phi[g, which_nonzero] * transient(pst, 
                                                                                               location = 0.75, scale = scale, reverse = reverse)
                                                   cells_on_constant_branch <- which(branch != which_nonzero)
                                                   cells_on_transient_branch <- which(branch == 
                                                                                        which_nonzero)
                                                   y <- rep(NA, C)
                                                   y[cells_on_transient_branch] <- rnorm(length(cells_on_transient_branch), 
                                                                                         mu[cells_on_transient_branch], gsd[g])
                                                   y[cells_on_constant_branch] <- X[g, cells_on_constant_branch]
                                                   return(y)
                                                 }))
  }
  if (zero_negative) {
    X[X < 0] <- 0
  }
  if (model_dropout && lambda < Inf) {
    drop_probs <- t(apply(X, 1, function(x) exp(-lambda * 
                                                  x)))
    for (g in seq_len(G)) {
      drop <- runif(C) < drop_probs[g, ]
      X[g, drop] <- 0
    }
  }
  X <- t(X)
  colnames(X) <- paste0("feature", 1:G)
  rownames(X) <- paste0("cell", 1:C)
  list(X = X, branch = branch, pst = pst, k = k, phi = phi, 
       delta = delta, p_transient = p_transient)
}


write.csv(synth[["X"]], "simulation.csv")
