PL.s <- function(theta, D.in, Z.in, obs.pairs = NULL, model = "MSP"){
  # theta (vector): dependence parameter, theta = log(c(sill, range))
  # D.in (matrix DxD): distance matrix for all sites
  # Z.in (matrix DxT): observations from corresponding sites
  # obs.pairs (matrix with 2 columns): selected observation (site) pairs
  # model (string): "MSP" or "IMSP"
  
  if (is.null(obs.pairs)){
    obs.pairs = t(do.call("cbind", sapply(1:(nrow(D.in)-1), function(k1){
      sapply((k1+1):nrow(D.in), function(k2){ c(k1,k2) } ) } ))) }
  
  npairs = nrow(obs.pairs)
  k1v = obs.pairs[,1]; k2v = obs.pairs[,2]
  
  D.pairs = sapply(1:nrow(obs.pairs), function(k) D.in[k1v[k], k2v[k]])
  a = sqrt(Vario.s(theta, D.pairs))

  if (model == "IMSP") {
    # Y.in = Z.in # IZ.in
    Z.in = 1/Z.in 
  } else if (model == "IMSP_trans") {
    Ztile.in = Z.in
    z1tile = Ztile.in[k1v,]; z2tile = Ztile.in[k2v,]
    Zstar.in = -log(1-exp(-1/Z.in))
    Z.in = 1/Zstar.in 
  }
  
  z1 = Z.in[k1v,]; z2 = Z.in[k2v,]
  log1 = log(z1) - log(z2); log2 = -log1
  
  u1 = a/2-log1/a; u2 = a/2-log2/a
  a = matrix(a, nrow = length(a), ncol = ncol(u1))
  
  # for mumerical stability
  trunc = 38
  case1 = which(u1 > trunc & u2 < -trunc)
  case2 = which(u1 < -trunc & u2 > trunc)
  case3 = which(u1 > trunc & u2 > trunc)
  case4ex = c(case1, case2, case3)
  
  if (model == "MSP") { 
    # max-stable processes
    if (length(case1) != 0) { 
      iz1 = 1/z1[case1]
      Cost1 = -sum(2*log(iz1) - iz1) 
    } else {Cost1 = 0}
    if (length(case2) != 0) { 
      iz2 = 1/z2[case2]
      Cost2 = -sum(2*log(iz2) - iz2) 
    } else {Cost2 = 0}
    if (length(case3) != 0) { 
      iz1 = 1/z1[case3]
      iz2 = 1/z2[case3]
      Cost3 = -sum(2*(log(iz1) + log(iz2)) - iz1 - iz2) 
    } else {Cost3 = 0}
    
    if (length(case4ex) != 0){
      a = a[-case4ex]
      u1 = u1[-case4ex]; u2 = u2[-case4ex]
      z1 = z1[-case4ex]; z2 = z2[-case4ex]
    }
    
    z1Sq = z1^2; z2Sq = z2^2
    az1z2 = a*z1*z2
    Phi1 = pnorm(u1); Phi2 = pnorm(u2)
    phi1 = dnorm(u1); phi2 = dnorm(u2)
    
    V = Phi1/z1 + Phi2/z2
    V1 = -Phi1/z1Sq - phi1/(a*z1Sq) + phi2/az1z2
    V2 = -Phi2/z2Sq - phi2/(a*z2Sq) + phi1/az1z2
    V12 = -(z2*u2*phi1 + z1*u1*phi2)/az1z2^2
    
    Cost4 = -sum(log(V1 * V2 - V12) - V)
    
  } else if (model %in% c("IMSP", "IMSP_trans")) { 
    # inverted max-stable processes
    Cost1 = Cost2 = 0
    
    if (length(case3) != 0) { 
      iz1 = 1/z1[case3]
      iz2 = 1/z2[case3]
      Cost3 = -sum(-(iz1 + iz2)) 
    } else {Cost3 = 0}
    
    if (length(case4ex) != 0){
      a = a[-case4ex]
      u1 = u1[-case4ex]; u2 = u2[-case4ex]
      z1 = z1[-case4ex]; z2 = z2[-case4ex]
      # a = a.rep[-case4ex]
    }
    
    z1Sq = z1^2; z2Sq = z2^2
    az1z2 = a*z1*z2
    Phi1 = pnorm(u1); Phi2 = pnorm(u2)
    phi1 = dnorm(u1); phi2 = dnorm(u2)
    
    V = Phi1/z1 + Phi2/z2
    V1 = -Phi1/z1Sq - phi1/(a*z1Sq) + phi2/az1z2
    V2 = -Phi2/z2Sq - phi2/(a*z2Sq) + phi1/az1z2
    V12 = -(z2*u2*phi1 + z1*u1*phi2)/az1z2^2
    
    Cost4 = -sum(log(z1Sq) + log(z2Sq) + log(V1 * V2 - V12) - V)
    
    if (model == "IMSP_trans") {
      b1 = -log(1-exp(-1/z1tile)) - 1/z1tile - 2*log(z1tile)
      b2 = -log(1-exp(-1/z2tile)) - 1/z2tile - 2*log(z2tile)
      Cost4 = Cost4 - sum(b1 + b2)
    }
  }
  
  Cost1 + Cost2 + Cost3 + Cost4
}


# Pairwise lig-likelihood for selected observation pairs (stationary model)
PL.ns <- function(theta, D.in, Z.in, reg.idx.in, pen.type = "None", 
                  adj.set = NULL, lambda = c(0, 0), obs.pairs = NULL,
                  model = "MSP"){
  # theta (vector): dependence parameter, theta = log(c(sill vector, range vector))
  # D.in (matrix DxD): distance matrix for all sites
  # Z.in (matrix DxT): observations from corresponding sites
  # reg.idx.in (vector D): partition indicies
  # pen.type (string): penalty used, L1, L2, or None
  # adj.set (matrix with 2 columns): neighboring subregion pairs
  # lambda (vector): penality parameters
  # obs.pairs (matrix with 2 columns): selected observation (site) pairs
  # model (string): "MSP" or "IMSP"
  
  if (is.null(obs.pairs)){
    obs.pairs = t(do.call("cbind", sapply(1:(nrow(D.in)-1), function(k1){
      sapply((k1+1):nrow(D.in), function(k2){ c(k1,k2) } ) } ))) }
  
  R.in = length(unique(as.vector(as.matrix(adj.set))))
  
  theta.mat = vec2mat(theta, R.in, lambda)
  
  sig.all = exp(theta.mat[1,])[reg.idx.in]
  phi.all = exp(theta.mat[2,])[reg.idx.in]
  
  k1v = obs.pairs[,1]; k2v = obs.pairs[,2]
  a.mat = sqrt(Vario.ns(sig.all, phi.all, D.in))
  a = sapply(1:nrow(obs.pairs), function(k) a.mat[k1v[k], k2v[k]])
  
  if (model == "IMSP") { 
    Z.in = 1/Z.in 
  } else if (model == "IMSP_trans") {
    Ztile.in = Z.in
    z1tile = Ztile.in[k1v,]; z2tile = Ztile.in[k2v,]
    Zstar.in = -log(1-exp(-1/Z.in))
    Z.in = 1/Zstar.in 
  }
  
  # ------------------------------------
  z1 = Z.in[k1v,]; z2 = Z.in[k2v,]
  log1 = log(z1) - log(z2); log2 = -log1
  
  u1 = a/2-log1/a; u2 = a/2-log2/a
  a = matrix(a, nrow = length(a), ncol = ncol(u1))
  
  trunc = 38
  case1 = which(u1 > trunc & u2 < -trunc)
  case2 = which(u1 < -trunc & u2 > trunc)
  case3 = which(u1 > trunc & u2 > trunc)
  case4ex = c(case1, case2, case3)
  
  if (model == "MSP") {
    if (length(case1) != 0) { 
      iz1 = 1/z1[case1]
      Cost1 = -sum(2*log(iz1) - iz1) 
    } else {Cost1 = 0}
    if (length(case2) != 0) { 
      iz2 = 1/z2[case2]
      Cost2 = -sum(2*log(iz2) - iz2) 
    } else {Cost2 = 0}
    if (length(case3) != 0) { 
      iz1 = 1/z1[case3]
      iz2 = 1/z2[case3]
      Cost3 = -sum(2*(log(iz1) + log(iz2)) - iz1 - iz2) 
    } else {Cost3 = 0}
    
    if (length(case4ex) != 0){
      a = a[-case4ex]
      u1 = u1[-case4ex]; u2 = u2[-case4ex]
      z1 = z1[-case4ex]; z2 = z2[-case4ex]
    }
    
    z1Sq = z1^2; z2Sq = z2^2
    az1z2 = a*z1*z2
    Phi1 = pnorm(u1); Phi2 = pnorm(u2)
    phi1 = dnorm(u1); phi2 = dnorm(u2)
    
    V = Phi1/z1 + Phi2/z2
    V1 = -Phi1/z1Sq - phi1/(a*z1Sq) + phi2/az1z2
    V2 = -Phi2/z2Sq - phi2/(a*z2Sq) + phi1/az1z2
    V12 = -(z2*u2*phi1 + z1*u1*phi2)/az1z2^2
    
    Cost4 = -sum(log(V1 * V2 - V12) - V)
  } else if (model %in% c("IMSP", "IMSP_trans")) {
    Cost1 = Cost2 = 0
    
    if (length(case3) != 0) { 
      iz1 = 1/z1[case3]
      iz2 = 1/z2[case3]
      Cost3 = -sum(-(iz1 + iz2)) 
    } else {Cost3 = 0}
    
    if (length(case4ex) != 0){
      a = a[-case4ex]
      u1 = u1[-case4ex]; u2 = u2[-case4ex]
      z1 = z1[-case4ex]; z2 = z2[-case4ex]
      # a = a.rep[-case4ex]
    }
    
    z1Sq = z1^2; z2Sq = z2^2
    az1z2 = a*z1*z2
    Phi1 = pnorm(u1); Phi2 = pnorm(u2)
    phi1 = dnorm(u1); phi2 = dnorm(u2)
    
    V = Phi1/z1 + Phi2/z2
    V1 = -Phi1/z1Sq - phi1/(a*z1Sq) + phi2/az1z2
    V2 = -Phi2/z2Sq - phi2/(a*z2Sq) + phi1/az1z2
    V12 = -(z2*u2*phi1 + z1*u1*phi2)/az1z2^2
    
    Cost4 = -sum(log(z1Sq) + log(z2Sq) + log(V1 * V2 - V12) - V)
    
    if (model == "IMSP_trans") {
      b1 = -log(1-exp(-1/z1tile)) - 1/z1tile - 2*log(z1tile)
      b2 = -log(1-exp(-1/z2tile)) - 1/z2tile - 2*log(z2tile)
      Cost4 = Cost4 - sum(b1 + b2)
    }
  }
  
  Cost = Cost1 + Cost2 + Cost3 + Cost4
  
  if (pen.type == 'L2'){
    pen = sum(sapply(1:length(lambda), function(i) {
      ifelse(lambda[i] == Inf, 0, lambda[i] * sum(sapply(1:nrow(adj.set), function(j){
        (theta.mat[i, adj.set[j,1]] - theta.mat[i, adj.set[j,2]])^2
      }))) }))
  } else if (pen.type == 'L1'){
    pen = sum(sapply(1:length(lambda), function(i) {
      ifelse(lambda[i] == Inf, 0, lambda[i] * sum(sapply(1:nrow(adj.set), function(j){
        abs(theta.mat[i, adj.set[j,1]] - theta.mat[i, adj.set[j,2]])
      }))) }))
  } else if (pen.type == 'None'){
    pen = 0
  }
  
  Cost + pen
}

