# Operate Partitions ###########################################################
minmaxStd = function(vec){ (vec-min(vec))/(max(vec)-min(vec)) }

Creat_par = function(R.in, S.in, type, seed_n = 1){
  # R.in: number of clusters/subregions
  # S.in (data frame with 2 columns): site coordinates
  stopifnot(is.integer(R.in))
  type = match.arg(type, c("grid", "cluster"))
  set.seed(seed_n)
  
  S.norm <- data.frame(lon = minmaxStd(S.in[,1]), lat = minmaxStd(S.in[,2]))
  
  # for grid
  if (R.in == 1L) {
    reg.idx = rep(1, nrow(S.in))
  } else {
    if (type == "grid") {
      coord = sapply(S.norm, function(x) {
        cuts <- seq(0,1,1/(sqrt(R.in)))
        cut(x, cuts, include.lowest = TRUE, labels = FALSE)
      })
      reg.idx = matrix(1:R.in, sqrt(R.in), sqrt(R.in))[coord]
    } else if (type == "cluster") {
      reg.idx = kmeans(S.norm, R.in, nstart = 25, iter.max = 10000)$cluster
    }
    
    reg.idx = sapply(reg.idx, function(id) which(unique(reg.idx) == id))
  }
  
  return(reg.idx)
}

# Find neighboring subregions
Neighbor = function(S, reg.idx){
  # S (data frame with 2 columns): site coordinates
  # reg.idx (vector D): subregion indicies representing the partition 
  
  S = as.data.frame(S)
  reg.uni = unique(reg.idx)
  R.in = length(reg.uni)
  A = matrix(0, nrow = R.in, ncol = R.in)
  colnames(A) = rownames(A) = reg.uni
  for (r in reg.uni){
    regs.idx = which(reg.idx == r)
    r.block = c()
    for (j in regs.idx){ # for each location
      p = as.numeric(S[j, ])
      for (k in 1:length(p)){
        p.col = S[S[,k]==p[k] & S[,3-k]!=p[3-k], ]
        nr = p.col[round(abs(p.col[,3-k] - p[3-k]),8)==round(min(abs(p.col[,3-k] - p[3-k])),8), ]
        r.block = c(r.block, reg.idx[as.numeric(rownames(nr))])
      }
    }
    r.block = unique(r.block)
    r.block = r.block[r.block!=r]
    r.block.pos = which(reg.uni %in% r.block)
    r.pos = which(reg.uni == r)
    A[r.pos, r.block.pos] = 1
  }
  A.upper = A; A.upper[lower.tri(A.upper)] = 0
  # adj.set = which(A.upper==1, arr.ind = T)
  adj.set = as.data.frame(which(A.upper==1, arr.ind = T))
  return(list(A = A, A.upper = A.upper, adj.set = adj.set))
}


# Creat Holdout & Validation, and Training Sets ################################
Creat_sets = function(S.in, Z.in, reg.idx, hold.perc, valid.perc, cv, seed_n, 
                      fold_n = 1L, ensure = T){
  
  set.seed(seed_n)
  
  if (cv) stopifnot(fold_n > 1L)
  stopifnot(is.integer(fold_n))
  
  numReg = length(unique(reg.idx))
  dataReg = lapply(1:numReg, function(r) {
    set.seed(seed_n)
    which(reg.idx==r)[shuffle(length(which(reg.idx==r)))]})
  allIdx = 1:length(reg.idx)
  
  if (!cv) { # single holdout set
    holdSize = floor(nrow(S.in)*hold.perc)
    
    # ensure at least one site in each subregion of base partition 
    if (ensure) {
      if (floor(holdSize/numReg) == 0) {
        holdIdx = sapply(sample(1:numReg, holdSize%%numReg), function(r) {dataReg[[r]][floor(holdSize/numReg)+1]})
      } else {
        holdIdx = c(as.vector(sapply(1:numReg, function(r) {dataReg[[r]][1:floor(holdSize/numReg)]})), 
                    sapply(sample(1:numReg, holdSize%%numReg), function(r) {dataReg[[r]][floor(holdSize/numReg)+1]}))
      }
    } else {
      holdIdx = sample((1:nrow(S.in)), holdSize)
    }
    
    validSize = floor(nrow(S.in)*valid.perc)
    validIdx = sample((1:nrow(S.in))[-holdIdx], validSize)
    
    # Validation set
    S.valid = S.in[validIdx, ]
    D.valid = fields::rdist(S.valid)
    Z.valid = Z.in[validIdx, ]
    
    # Training set for stationary model
    trainIdx.s = setdiff(allIdx, validIdx)#allIdx[-validIdx]
    S.train.s = S.in[trainIdx.s, ]     # training sites for stationary model
    D.train.s = fields::rdist(S.train.s)       # corresponding distance matrix
    Z.train.s = Z.in[trainIdx.s, ]                # corresponding data
    
    # Training and holdout sets for nonstationary model
    S.hold = S.in[holdIdx, ] 
    D.hold = fields::rdist(S.hold)
    Z.hold = Z.in[holdIdx, ]
    
    trainIdx.ns = setdiff(allIdx, c(holdIdx, validIdx))
    S.train.ns = S.in[trainIdx.ns, ]
    D.train.ns = fields::rdist(S.train.ns)
    Z.train.ns = Z.in[trainIdx.ns, ]
    
    
  } else if (cv) {
    # data fold
    # holdout one site for each subregion
    dataReg.size = sapply(1:length(dataReg), function(r) { length(dataReg[[r]]) })
    
    if (ensure & min(dataReg.size) >= fold_n) {
      validSize = floor(nrow(S.in)*valid.perc)
      
      holdIdx1 = lapply(1:fold_n, function(i) {
        sapply(1:numReg, function(r){
          dataReg[[r]][i]
        })})
      restId = unlist(sapply(1:numReg, function(r){
        dataReg[[r]][(fold_n+1):length(dataReg[[r]])]
      }))
      
      # validation set
      validIdx = sample(restId, validSize)
      
      # exclude validation set
      restId = restId[which(!(restId %in% validIdx))]
      holdIdx2i = caret::createFolds(restId, k=fold_n, list=T)
      
      # holdout set
      # each fold contains at least one location in each subregion in base partition
      holdIdx = lapply(1:fold_n, function(i){
        c(holdIdx1[[i]], restId[holdIdx2i[[i]]])
      })
    } else {
      validSize = floor(nrow(S.in)*valid.perc)
      validIdx = sample(1:nrow(S.in), validSize)
      holdSize = floor(nrow(S.in)*hold.perc)
      restId = (1:nrow(S.in))[-validIdx]
      cut_fold = cut(shuffle(restId),breaks=fold_n,labels=FALSE)
      # restId[shuffle(restId)]
      holdIdx = lapply(1:fold_n, function(i){
        restId[which(cut_fold == i)]
      })
    }
    
    # Validation
    S.valid = S.in[validIdx, ]
    D.valid = fields::rdist(S.valid)
    Z.valid = Z.in[validIdx, ]
    
    # Training set for stationary model
    trainIdx.s = allIdx[-validIdx]
    
    S.train.s = S.in[trainIdx.s, ]     # training sites for stationary model
    D.train.s = fields::rdist(S.train.s)       # corresponding distance matrix
    Z.train.s = Z.in[trainIdx.s, ]                # corresponding data
    
    # Training and holdout sets for nonstationary model
    S.hold = lapply(1:fold_n, function(i) S.in[holdIdx[[i]], ] )
    D.hold = lapply(1:fold_n, function(i) fields::rdist(S.hold[[i]]) )
    Z.hold = lapply(1:fold_n, function(i) Z.in[holdIdx[[i]], ] )
    
    trainIdx.ns = lapply(1:fold_n, function(i) setdiff(allIdx, c(holdIdx[[i]], validIdx)) )
    S.train.ns = lapply(1:fold_n, function(i) S.in[trainIdx.ns[[i]], ] )
    D.train.ns = lapply(1:fold_n, function(i) fields::rdist(S.train.ns[[i]]) )
    Z.train.ns = lapply(1:fold_n, function(i) Z.in[trainIdx.ns[[i]], ] )
    
  }
  
  list(validIdx = validIdx,
       S.valid = S.valid, 
       D.valid = D.valid, 
       Z.valid = Z.valid,
       trainIdx.s = trainIdx.s, 
       S.train.s = S.train.s,
       D.train.s = D.train.s,
       Z.train.s = Z.train.s, 
       holdIdx = holdIdx,
       S.hold = S.hold,
       D.hold = D.hold,
       Z.hold = Z.hold,
       trainIdx.ns = trainIdx.ns,
       S.train.ns = S.train.ns,
       D.train.ns = D.train.ns,
       Z.train.ns = Z.train.ns)
}


# Sample Observation Pairs #####################################################
Samp = function(scheme, obs.pairs, frac, dist_mat, num.itv = 1L) {
  # scheme (string): "sr"-simple random sampling, "st"-stratified samplings
  # obs.pairs (matrix with 2 columns): all observation pairs
  # dist_mat (matrix DxD): distance matrix for all sites 
  # num.itv (integer): number of stratification
  stopifnot(is.integer(num.itv))
  scheme = match.arg(scheme, c("sr", "st"))
  
  if (scheme == "sr") {
    sample.pairs = sample(1:nrow(obs.pairs), round(nrow(obs.pairs)*frac))
  } else if (scheme == "st") {
    stopifnot(num.itv >= 1L & is.integer(num.itv))
    
    # Distance for observation pairs
    dist.pairs = sapply(1:nrow(obs.pairs), function(r){
      pair = obs.pairs[r,]; dist_mat[pair[1], pair[2]] })
    divide = seq(min(dist.pairs),max(dist.pairs),length.out = num.itv + 1)
    divide[length(divide)] = divide[length(divide)] + 0.1
    divide.idx = lapply(1:num.itv, function(i) {
      which(divide[i] <= dist.pairs & dist.pairs < divide[i+1])
    })
    sample.pairs = c()
    for (i in 1:num.itv) {
      sample.pairs = c(sample.pairs, sample(divide.idx[[i]], max(1,round(length(divide.idx[[i]])*pair.frac))))
    }
  }
  
  sample.pairs
}

# Select observation pairs
SelPairs = function(idx.vec, frac, dist_mat, seed_n, type = "st", itv_n = 20L) {
  set.seed(seed_n)
  obsPairs = t(do.call("cbind", sapply(1:(length(idx.vec)-1), function(k1){
    sapply((k1+1):length(idx.vec), function(k2){ c(k1,k2) } ) } )))
  samplePairs = Samp(type, obsPairs, frac, dist_mat, itv_n)
  obsPairs[samplePairs, ]
}


# Pairwise Log-Likelihood ######################################################
# # Creat Observation Pairs
# CreatPairs = function(idx.vec) {
#   t(do.call("cbind", sapply(1:(length(idx.vec)-1), function(k1){
#     sapply((k1+1):length(idx.vec), function(k2){ c(k1,k2) } ) } )))
# }

# Stationary variogram
Vario.s = function(theta, D.in){ 
  # theta (vector): dependence parameter, theta = c(sill, range)
  # D.in (matrix): distance matrix of sites
  2*exp(theta[1]) - 2*exp(theta[1])*exp(-D.in/exp(theta[2]))
}

# Nonstationary variogram
Vario.ns = function(sig.all, phi.all, D.in){
  # sig.all (vector): a vector containing sill estimates for all sites
  # phi.all (vector): a vector containing range estimates for all sites
  # D.in (matrix DxD): distance matrix for all sites
  n = nrow(D.in)
  sig.sum = matrix(sig.all, n, n) + matrix(sig.all, n, n, byrow = T)
  sig.mul = sqrt(sig.all)%*%t(sqrt(sig.all))
  Rns1 = phi.all%*%t(phi.all)
  Rns2 = 2 / (matrix(phi.all^2, n, n) + 
                matrix(phi.all^2, n, n, byrow = T))
  Rns = Rns1*Rns2*exp(-sqrt(Rns2)*D.in)
  double.gamma = sig.sum - 2 * sig.mul * Rns
  diag(double.gamma) = 0
  return(double.gamma)
}

# RI, LRI ######################################################################

RI = function(P0, P1){
  if (length(P0) == length(P1)) {
    D = length(P0)
  } else { return(0) }
  dom = D*(D-1)/2
  num_ss = 0; num_dd = 0
  for (i in 1:(D-1)){
    for (j in (i+1):D) {
      if (P0[i] == P0[j] & P1[i] == P1[j]) {
        num_ss = num_ss + 1
      } else if (P0[i] != P0[j] & P1[i] != P1[j]) {
        num_dd = num_dd + 1
      }
    }
  }
  (num_ss + num_dd) / dom
}


LRI = function(P0, P1){
  # for each partition P1, true partition P0
  if (length(P0) == length(P1)) {
    D = length(P0)
  } else { return(0) }
  beta.all = c()
  for (i in 1:D){
    betai = 0
    num_ss = 0; num_dd = 0
    for (j in (1:D)[-i]) {
      if (P0[i] == P0[j] & P1[i] == P1[j]) {
        num_ss = num_ss + 1
      } else if (P0[i] != P0[j] & P1[i] != P1[j]) {
        num_dd = num_dd + 1
      }
    }
    betai = num_ss + num_dd
    beta.all = c(beta.all, betai)
  }
  return(beta.all/D)
}


# Others #######################################################################
Get_inival = function(est.old, reg.idx.new, reg.idx.old, lambda = c(0, 0)){
  R.old = length(unique(reg.idx.old))
  R.new = length(unique(reg.idx.new))
  ini.val = c()
  idx = 1
  # ---
  for (k in 1:length(lambda)){
    if (lambda[k] == Inf){
      # ini.val = c(ini.val, est.old[idx])
      # idx = idx + 1
      ini.val = c(ini.val, mean(est.old[idx:(idx+R.old-1)]))
      idx = idx + R.old
    } else {
      est.all = est.old[idx:(idx+R.old-1)][reg.idx.old]
      extract = sapply(1:R.new, function(r){mean(est.all[reg.idx.new==r])})
      ini.val = c(ini.val, extract)
      idx = idx + R.old
    }
  }
  # ---
  ini.val
}

trans = function(theta, R.in = 1L, lambda = c(0,0), repara = c("log", "log"), 
                 bounds = c(2, 2)){
  lens = c(sum(lambda[1] == Inf) + R.in*sum(lambda[1] != Inf), 
           sum(lambda[2] == Inf) + R.in*sum(lambda[2] != Inf))
  theta.list = list(theta[1:lens[1]], theta[lens[1]+1:lens[2]])
  theta.trans = c()
  for (k in 1:length(repara)){
    if (repara[k] == "log"){
      theta.trans = c(theta.trans, log(theta.list[[k]]))
    } 
    else if (repara[k] == "logit") {
      theta.trans = c(theta.trans, logit(theta.list[[k]]/bounds[k]))
    }
  }
  return(theta.trans)
}

invtrans = function(theta, R.in = 1L, repara = c("log", "log"), 
                    lambda = c(0,0), bounds = c(2, 2)){
  lens = c(sum(lambda[1] == Inf) + R.in*sum(lambda[1] != Inf), 
           sum(lambda[2] == Inf) + R.in*sum(lambda[2] != Inf))
  theta.list = list(theta[1:lens[1]], theta[lens[1]+1:lens[2]])
  theta.trans = c()
  for (k in 1:length(repara)){
    if (repara[k] == "log"){
      theta.trans = c(theta.trans, exp(theta.list[[k]]))
    } 
    else if (repara[k] == "logit") {
      theta.trans = c(theta.trans, bounds[k]*invlogit(theta.list[[k]]))
    }
  }
  return(theta.trans)
}

vec2mat = function(vec, R.in, lambda) {
  mat = matrix(0,length(lambda),R.in); idx = 1
  for (k in 1:length(lambda)){
    if (lambda[k] == Inf){
      mat[k,] = rep(vec[idx],R.in); idx = idx + 1
    } else { mat[k,] = vec[idx:(idx+R.in-1)]; idx = idx + R.in }}
  mat
}

# mat2vec = function() {}

autoStopCluster <- function(cl) {
  # Acknowledgement: https://stackoverflow.com/questions/52190651/how-to-shut-down-an-open-r-cluster-connection-using-parallel
  stopifnot(inherits(cl, "cluster"))
  env <- new.env()
  env$cluster <- cl
  attr(cl, "gcMe") <- env
  reg.finalizer(env, function(e) {
    message("Finalizing cluster ...")
    message(capture.output(print(e$cluster)))
    try(parallel::stopCluster(e$cluster), silent = FALSE)
    message("Finalizing cluster ... done")
  })
  cl
}

