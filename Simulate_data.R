library(mvtnorm)
library(permute)
library(fields)
library(LaplacesDemon)
library(plot3D)
library(ggplot2)
library(parallel)
library(doParallel)
require(this.path)

# RNGkind(sample.kind = "Rejection")
setwd(this.path::here())

source('Functions/Utils.R')
seed_n = 1
set.seed(seed_n)

################################################################
# Data generation ##############################################
################################################################

ngrid = 20
grids = seq(0, 1, length.out = ngrid)
S <- expand.grid(x1 = grids, x2 = grids)    # location grid
npro <- nrow(S)                             # number of locations
D = rdist(S)                                # distance matrix
B <- 100000                                 # block size
Trep = 100                                  # number of replications

# A Simple 2x2 square grid #####################################

partition = sapply(S[,c('x1','x2')], function(x) {
  cuts <- c(0, 0.5, 1)
  cut(x, cuts, include.lowest = TRUE, labels = FALSE)
})

# Varying Sill & Fixed Range ###################################

sig.vec = c(0.5, 2, 2, 8)
sig.all = matrix(sig.vec, 2, 2, byrow = T)[partition]
phi = 0.2; phi.vec = rep(phi, 4)
phi.all = matrix(phi.vec, 2, 2, byrow = T)[partition]

# Varying Range & Fixed Sill ###################################

# sigSq = 2; sig.vec = rep(sigSq, 4)
# sig.all = matrix(sig.vec, 2, 2, byrow = T)[partition]
# phi.vec = c(0.05, 0.2, 0.2, 0.8)
# phi.all = matrix(phi.vec, 2, 2, byrow = T)[partition]

# Nonstationary Variogram  #####################################

vario.ns = Vario.ns(sig.all, phi.all, D)
covmat.ns = (matrix(sig.all, npro, npro) + 
               matrix(sig.all, npro, npro, byrow = T))/2 - vario.ns/2
halfVar = sig.all/2

# Data Generation  #############################################

# Speed up with parallel computation 
ncores = 5
cl <- makeCluster(ncores)
registerDoParallel(cl)
clusterExport(cl, varlist = c("Trep", "B", "covmat.ns", "halfVar", "npro"))
invisible(clusterEvalQ (cl , library("mvtnorm"))) 

Z.max <- parSapply(cl, 1:Trep, function(n){
  PPPn <- cumsum(rexp(B))
  iPPPn = -log(PPPn)
  GPn = t(rmvnorm(n=B, mean=rep(0, npro), sigma=covmat.ns, method="chol"))
  logZn = apply(sapply(1:B, function(b){GPn[,b] - halfVar + iPPPn[b]}), 1, max)
  return(exp(logZn))})

stopCluster(cl)

save(S, Z.max, partition, sig.all, phi.all, file=paste0("Simulated.Rdata"))
