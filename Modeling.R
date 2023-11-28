require(mvtnorm)
require(permute)
require(fields)
require(LaplacesDemon)
require(ggplot2)
require(parallel)
require(doParallel)
require(caret)
require(viridis)
require(this.path)

rm(list = ls())

################################################################################
# Basic Setting ################################################################
################################################################################

# Input Data ###################################################################
# input renormalized data
setwd(this.path::here())
dat = "Simulated" # "NepalExtended" # 
load(paste0(dat, ".Rdata"))

model = "MSP"


if (model == "IMSP") {
  iZ.max = 1/Z.max
  Z.max = iZ.max
}
S = data.frame(S)           # site coordinates
str(Z.max)                  # renormalized data

source('Functions/Utils.R')           
source('Functions/Fit.R')
source('Functions/Objectives.R')
source('Functions/Lambda_tuning.R')
source('Functions/Merge_subr.R')
source('Functions/Algorithm1.R')

seed = 1                    # seed number for sampling
# Hyper-parameter ##############################################################
if (dat == "Simulated"){
  # For simulated data
  # Modeling of this simulated data takes around ten minutes
  pen.type = "L1"           # "L1" or "L2"
  BP.type = "grid"          # the way we generate base partition: "cluster" or "grid" (the latter only works for the square domain)
  ncores = 3                # number of cores for parallel computation
  cv = F                    # "T" for cross-validation-type or "F" for single holdout set
  fold_n = 5L               # number of data fold (for data fold case)
  pair.frac = 1/100         # fraction of observation pairs. Smaller dataset should consider larger fraction.
  R.base = 64L              # number of subregions in base partition
  Lambda = c(Inf, 2^c(5,4,3,2,1,0,-1), 0) # descending grid, may change according to data
  stp.len = 0               # step length for the partition merging search, 0 represents the self-adapted step length
  layers = c(4L,16L)        # hierarchical structure for acquiring proper starting point
} else if (dat == "NepalExtended"){
  # For extended Nepal data
  # !!!WARNING: The modeling for the extended Nepal dataset might take 4-6 hours in total
  pen.type = "L2"           # "L1" or "L2"
  BP.type = "cluster"       # the way we generate base partition: "cluster" or "grid" (the latter only works for the square domain)
  ncores = 5                # number of cores for parallel computation
  cv = T                    # "T" for cross-validation-type or "F" for single holdout set
  fold_n = 5L               # number of data fold (for data fold case)
  pair.frac = 1/500         # fraction of observation pairs
  R.base = 80L              # number of subregions in base partition
  Lambda = c(Inf, 2^c(5,4,3,2,1,0,-1), 0) # descending grid, may change according to data
  stp.len = 0.2             # step length for the partition merging search, 0 represents the self-adapted step length
  layers = c(8L,20L,40L)    # hierarchical structure for acquiring proper starting point
}

RNGkind(sample.kind = "Rejection")

set.seed(seed)
# Creat Base Partition #########################################################
# generate base subregion partition with given number
reg.idx.base = Creat_par(R.base, S, BP.type, seed)
plot(S, col = sample(colors(), R.base)[reg.idx.base], pch = 15)

# Holdout & Validation, and Training Sets ######################################
sets = Creat_sets(S, Z.max, reg.idx.base, hold.perc = 0.15, valid.perc = 0.15, 
                  cv, seed_n = seed, fold_n = fold_n)

################################################################################
# Stationary Model #############################################################
################################################################################
# Sampled observation pairs for stationary model
samplePairs.s = SelPairs(sets$trainIdx.s, pair.frac, sets$D.train.s,
                         itv_n = 20L, seed_n = seed)

mod.s0 = Fit.s(c(0.5, 0.5), sets$D.train.s, sets$Z.train.s, samplePairs.s, model = model)
est.s0 = mod.s0$est
pl.s0 = -PL.s(mod.s0$est.trans, sets$D.valid, sets$Z.valid, model = model)

cat(" Parameter estimates with stationary model:", round(est.s0[1],3), 
    ' (sill) and', round(est.s0[2], 3), ' (range).\n', 
    "Optimization time:", as.numeric(mod.s0$time), units(mod.s0$time), "\n",
    "Pairwise likelihood for stationary model:", round(pl.s0))

model.s = list(est.s = est.s0, pl.s = pl.s0)

################################################################################
# Pre-estimation ###############################################################
################################################################################
# Estimate a proper initial estimates for nonstationary models, 
# which is essential for high dimensional data
print("Initial Estimates")

if (cv) {
  cl <- autoStopCluster(makeCluster(ncores))
  # cl <- makeCluster(ncores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("R.base", "S", "reg.idx.base", "reg.idx.base",
                                "pair.frac", "sets", "pen.type", 
                                "cv", "fold_n", "model",
                                "Fit.s", "Fit.ns", "PL.s", "PL.ns",
                                "Vario.ns", "Vario.s", "trans", "invtrans",
                                "Creat_par", "minmaxStd", "Get_inival",
                                "Neighbor", "vec2mat", "rdist",
                                "SelPairs", "Samp"))
}


ini.est_file = paste0("Modeling/", dat, "_", model, "_ini_seed", seed, "_", cv, 
                      "_", BP.type, "_", paste(layers, collapse = "-"), ".Rdata")

if (file.exists(ini.est_file)) {
  load(ini.est_file)
  print("Input existing StartPtFile, which includes a proper starting point for subsequent optimization.")
} else { 
  # set.seed(seed)
  model.ini = Fit.ini(layers = layers, S, 
                      sets, cv, pair.frac, fold_n, 
                      BP.type, seed = seed, model = model) 
  save(model.ini, file = ini.est_file)
}

print(model.ini$time)

est.s = model.ini$est.s
if (!cv) {pl.s = model.ini$pl.s} else if (cv) {pl.s = model.ini$pl.s.mean}
samplePairs.ns = model.ini$samplePairs.ns

para = model.ini$est.ini[1:model.ini$R.ini][model.ini$reg.idx.ini]
ggplot(data.frame(S), aes(x = data.frame(S)[ , 1], y = data.frame(S)[ , 2], fill = para)) +
  geom_tile() + scale_fill_viridis() + xlab("S1") +  ylab("S2") +
  ggtitle("Initial Values (Sill)")
para = model.ini$est.ini[1:model.ini$R.ini+model.ini$R.ini][model.ini$reg.idx.ini]
ggplot(data.frame(S), aes(x = data.frame(S)[ , 1], y = data.frame(S)[ , 2], fill = para)) +
  geom_tile() + scale_fill_viridis() + xlab("S1") +  ylab("S2") +
  ggtitle("Initial Values (Range)")




################################################################################
# Fit Nonstationary Model for Base Partition ###################################
################################################################################
Lambda.list = list(Lambda,Lambda)
lambda = c(Inf, Inf) # initial lambda corresponding stationary model

# Lambda Tunning ###############################################################
model.base = lamTuning(lambda, Lambda.list, loc = S,
                       reg.idx.in = reg.idx.base, est.in = est.s, pl.in = pl.s, 
                       pen.type = pen.type, sets = sets, mod.ini = model.ini, 
                       cv = cv, show = T, ncores = ncores, model = model)
# Input lambda
print(model.base$lambda.in) 
# Optimal lambda for base partition
print(model.base$lambda) 
# Lambda-tuning time (min)
print(model.base$time/60) 


################################################################################
# Merge subregions #############################################################
################################################################################

model.merge = Alg1(R.in = R.base, 
                   reg.idx.in = reg.idx.base, 
                   mod.ini = model.ini, 
                   mod.base = model.base, 
                   loc = S,
                   sets = sets, 
                   pen.type = pen.type, 
                   fold_n = fold_n, 
                   stp.len = stp.len, 
                   model = model)

# Subregion-merging time (min)
print(model.merge$RT/60)


# Save Model ###################################################################
models.save = list(hyper = list(pen.type = pen.type,
                                BP.type = BP.type,
                                ncores = ncores,
                                seed = seed,
                                cv = cv,
                                fold_n = fold_n,
                                pair.frac = pair.frac,
                                R.base = R.base,
                                Lambda = Lambda,
                                layers = layers,
                                stp.len = stp.len),
                   sets = sets, model.s = model.s, 
                   model.ini = model.ini, model.base = model.base, 
                   model.merge = model.merge)
save(models.save, file = paste0("Modeling/", dat, "_", model, "_", pen.type, 
                                "_seed", seed, "_", cv, "_Results.Rdata"))


