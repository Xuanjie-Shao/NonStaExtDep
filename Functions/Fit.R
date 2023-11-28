# Some fitting functions for the convenience
# Estimation ###################################################################
Fit.s = function(par, D.in, Z.in, obs.pairs, hessian = F, model ="MSP") {
  ### Fit the stationary model
  # par (vector): dependence parameter, par = c(sill, range)
  # D.in (matrix DxD): distance matrix for all sites
  # Z.in (matrix DxT): observations from corresponding sites
  # obs.pairs (matrix with 2 columns): selected observation (site) pairs
  # model (string): "MSP" or "IMSP"
  tmp = Sys.time()
  mod = optim(par = trans(par), PL.s, hessian = hessian,
              D.in = D.in, Z.in = Z.in, obs.pairs = obs.pairs,
              model = model)
  est = invtrans(mod$par)
  
  list(mod = mod, est = est, est.trans = mod$par, time = Sys.time() - tmp)
}

Fit.ns = function(par, D.in, Z.in, R.in, reg.idx.in, adj.set, pen.type, obs.pairs, 
                  hessian = F, lambda = c(0,0), model ="MSP") {
  ### Fit the nonstationary model
  # par (vector): dependence parameter, par = c(sill vector, range vector)
  # D.in (matrix DxD): distance matrix for all sites
  # Z.in (matrix DxT): observations from corresponding sites
  # reg.idx.in (vector D): partition indicies
  # adj.set (matrix with 2 columns): neighboring subregion pairs
  # pen.type (string): penalty used, L1, L2, or None
  # obs.pairs (matrix with 2 columns): selected observation (site) pairs
  # lambda (vector): penality parameters
  # model (string): "MSP" or "IMSP"
  
  tmp = Sys.time()
  mod = optim(par = trans(par, R.in, lambda = lambda), PL.ns, hessian = hessian,
             D.in = D.in, Z.in = Z.in, reg.idx.in = reg.idx.in, adj.set = adj.set,
             pen.type = pen.type, obs.pairs = obs.pairs, lambda = lambda,
             model = model)
  est = invtrans(mod$par, R.in, lambda = lambda)
  
  list(mod = mod, est = est, est.trans = mod$par, time = Sys.time() - tmp)
}

# Pre-Estimation ###############################################################
Fit.ini = function(layers = NULL, loc, sets, cv, pair.frac, fold_n, BP.type, 
                   pen.type = "None", seed = 1, model ="MSP", est.s = NULL) {
  ### Fit the initial valuess
  tmp = Sys.time()
  if (!cv) {

    if (is.null(est.s)) {
      samplePairs.ns = SelPairs(sets$trainIdx.ns, pair.frac, sets$D.train.ns, seed)
      mod.s = Fit.s(c(0.5, 0.5), sets$D.train.ns, sets$Z.train.ns, samplePairs.ns, 
                    model = model)
      est.s = mod.s$est
      est.s.trans = mod.s$est.trans
    } else {
      est.s.trans = trans(est.s)
    }
    pl.s = -PL.s(est.s.trans, sets$D.hold, sets$Z.hold, model = model)
    
    R.ini = 1L
    reg.idx.ini = rep(1, nrow(loc))
    est.ini = est.s
    if (!is.null(layers)) {
      for (j in 1:length(layers)){
        cat("Layers:", layers[j], "\n")
        R0 = R.ini
        reg.idx0 = reg.idx.ini
        est0 = est.ini
        
        R.ini = layers[j]
        reg.idx.ini = Creat_par(R.ini, loc, BP.type, seed_n = seed)
        ini.val = Get_inival(est0, reg.idx.new = reg.idx.ini, reg.idx.old = reg.idx0)
        
        adj.set.ini = Neighbor(loc, reg.idx.ini)$adj.set
        mod.ns.ini = Fit.ns(ini.val, sets$D.train.ns, sets$Z.train.ns, 
                            R.ini, reg.idx.ini[sets$trainIdx.ns],
                            adj.set.ini, pen.type, samplePairs.ns, hessian = T, 
                            model = model)
        est.ini = mod.ns.ini$est
      }
    }
    
    list(samplePairs.ns = samplePairs.ns, est.s = est.s, pl.s = pl.s,
         est.ini = est.ini, R.ini = R.ini, reg.idx.ini = reg.idx.ini, 
         time = Sys.time() - tmp)
  } else if (cv) {
    # Fit stationary model for each nonstationary training set
    Model = foreach(i = 1:fold_n) %dopar% {
      # Selected data fold will be set as a holdout set
      # The rest will be combined as a nonstationary training set

      if (is.null(est.s)) {
        samplePairs.ns = SelPairs(sets$trainIdx.ns[[i]], pair.frac, sets$D.train.ns[[i]], seed)
        mod.s = Fit.s(c(0.5, 0.5), sets$D.train.ns[[i]], sets$Z.train.ns[[i]], samplePairs.ns, 
                      model = model)
        est.s = mod.s$est
        est.s.trans = mod.s$est.trans
      } else {
        est.s.trans = trans(est.s)
      }
      
      pl.s = -PL.s(mod.s$est.trans, sets$D.hold[[i]], sets$Z.hold[[i]], model = model)
      
      # -----
      R.ini = 1L
      reg.idx.ini = rep(1, nrow(loc))
      est.ini = est.s
      if (!is.null(layers)) {
        for (j in 1:length(layers)){
          cat("Layers:", layers[j])
          R0 = R.ini
          reg.idx0 = reg.idx.ini
          est0 = est.ini
          
          R.ini = layers[j]
          reg.idx.ini = Creat_par(R.ini, loc, BP.type, seed_n = seed)
          ini.val = Get_inival(est0, reg.idx.new = reg.idx.ini, reg.idx.old = reg.idx0)
          
          adj.set.ini = Neighbor(loc, reg.idx.ini)$adj.set
          mod.ns.ini = Fit.ns(ini.val, sets$D.train.ns[[i]], sets$Z.train.ns[[i]], 
                              R.ini, reg.idx.ini[sets$trainIdx.ns[[i]]],
                              adj.set.ini, pen.type, samplePairs.ns, hessian = T,
                              model = model)
          est.ini = mod.ns.ini$est
        }
      }
      # -----
      
      list(samplePairs.ns = samplePairs.ns,
           est.s = est.s, pl.s = pl.s, 
           est.ini = est.ini) 
    }

    samplePairs.ns = lapply(1:fold_n, function(i) Model[[i]]$samplePairs.ns)
    est.s = lapply(1:fold_n, function(i) {Model[[i]]$est.s})
    pl.s = lapply(1:fold_n, function(i) Model[[i]]$pl.s)
    est.s.mean = rowMeans(sapply(1:fold_n, function(i) {Model[[i]]$est.s}))
    pl.s.mean = mean(sapply(1:fold_n, function(i) Model[[i]]$pl.s))
    
    R.ini = ifelse(is.null(layers), 1L, layers[length(layers)])
    est.ini = rowMeans(sapply(1:fold_n, function(i) {Model[[i]]$est.ini}))
    reg.idx.ini = Creat_par(R.ini, loc, BP.type, seed_n = seed)
    
    list(samplePairs.ns = samplePairs.ns, est.s = est.s, pl.s = pl.s,
         est.s.mean = est.s.mean, pl.s.mean = pl.s.mean,
         est.ini = est.ini, R.ini = R.ini, reg.idx.ini = reg.idx.ini, 
         time = Sys.time() - tmp)
  }
  
}



