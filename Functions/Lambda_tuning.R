lamTuning = function(lambda, Lambda.list, loc, reg.idx.in, est.in, pl.in, 
                     pen.type, sets, mod.ini, cv = FALSE, show = FALSE, 
                     ncores = 5L, model = "MSP"){
  
  # lambda (vector): penalty parameter
  # Lambda.list (list with 2 vectors): descending grid for lambda-tuning
  # loc (matrix D*2): coordinate of the sites
  # reg.idx.in (vector D): partition with subregion indicies
  # est.in (vector): parameter estimates with current lambda
  # pl.in (float): pairwise likelihood with current parameter estimates
  # sets (list): the setting of traning set, holdout set, validation set, and so on
  # pen.type (string): penalty used, "L1", "L2", or "None"
  # mod.ini (list): modeling results for initial values
  # cv (boolean): cross-validation-type holdout sets or single holdout set
  # show (logistic): show iteration information or not
  # ncores (integer): number of the cores for parallel computing
  # model (string): "MSP" or "IMSP"
  
  cat("----- Lambda Tuning ----- \n")
  
  ptm = proc.time()
  
  est.ini = mod.ini$est.ini
  reg.idx.ini = mod.ini$reg.idx.ini
  obs.pairs = mod.ini$samplePairs.ns
  
  R.in = length(unique(reg.idx.in))
  # A.upper = Neighbor(loc, reg.idx)$A.upper
  # adj.set.pos = which(A.upper==1, arr.ind = T)
  adj.set.in = Neighbor(loc, reg.idx.in)$adj.set
  # est.trans = trans(est.in, R.in, lambda)
  lambda.in = lambda
  
  est.tune = est.in
  pl.tune = pl.in
  
  
  
  tuning = TRUE
  while (tuning) {
    pl.collect = c() #
    lambda.collect = matrix(NaN, nrow = length(lambda), ncol = length(lambda)) #
    eta.collect = list() #
    lambda.former = lambda ##
    for (k in 1:length(lambda)){
      lambda = lambda.former
      if (lambda[k] != Lambda.list[[k]][length(Lambda.list[[k]])]){
        # uodate lambda
        lambda[k] = Lambda.list[[k]][which(Lambda.list[[k]] == lambda[k])+1]
        lambda.collect[k, ] = lambda
        
        ##### ---------------------------
        if (!cv) {
          ini.val = Get_inival(est.ini, reg.idx.in, reg.idx.ini, lambda)
          mod.tmp = Fit.ns(ini.val, sets$D.train.ns, sets$Z.train.ns, 
                           R.in, reg.idx.in[sets$trainIdx.ns],
                           adj.set.in, pen.type, obs.pairs, 
                           lambda = lambda, model = model) 
          
          est.tmp = mod.tmp$est
          # eta.collect = c(eta.collect, list(est.tmp))
          pl.tmp = -PL.ns(mod.tmp$est.trans, sets$D.hold, sets$Z.hold, reg.idx.in[sets$holdIdx], 
                        pen.type = 'None', adj.set = adj.set.in, lambda = lambda, model = model)
        } else if (cv) {
          Model.tmp = foreach(i = 1:fold_n) %dopar% {
            ini.val = Get_inival(est.ini, reg.idx.in, reg.idx.ini, lambda) #[[i]]
            mod.tmp = Fit.ns(ini.val, sets$D.train.ns[[i]], sets$Z.train.ns[[i]], 
                             R.in, reg.idx.in[sets$trainIdx.ns[[i]]],
                             adj.set.in, pen.type, obs.pairs[[i]], 
                             lambda = lambda, model = model)
            
            est.tmp = mod.tmp$est
            # eta.collect = c(eta.collect, list(est.tmp))
            pl.tmp = -PL.ns(mod.tmp$est.trans, sets$D.hold[[i]], sets$Z.hold[[i]], 
                          reg.idx.in[sets$holdIdx[[i]]], pen.type = 'None', 
                          adj.set = adj.set.in, lambda = lambda, model = model)
            list(est.tmp = est.tmp, pl.tmp = pl.tmp)
          }
          
          est.tmp = lapply(1:fold_n, function(i) Model.tmp[[i]]$est.tmp)
          # eta.collect = c(eta.collect, list(est.tmp))
          pl.tmp = mean(sapply(1:fold_n, function(i) Model.tmp[[i]]$pl.tmp))
        }
        ##### ---------------------------
        
        eta.collect = c(eta.collect, list(est.tmp))
        pl.collect = c(pl.collect, pl.tmp) #
      } else if (lambda[k] == Lambda.list[[k]][length(Lambda.list[[k]])]){
        lambda.collect[k, ] = lambda.former
        eta.collect = c(eta.collect, list(est.tune))
        pl.collect = c(pl.collect, pl.tune)
      }
    }
    
    if (max(pl.collect) - pl.tune > 0){
      max.id = which(pl.collect == max(pl.collect))[1]
      lambda = lambda.collect[max.id, ] ###
      est.tune = eta.collect[[max.id]] ###
      if (show == T){
        print("# -----")
        cat("Tuned lambda:", paste0(lambda, collapse = ", "), "\n")
        cat("Current pairwise log-likelihood:", pl.tune, "\n")
        cat("Maximum tried pairwise log-likelihood:", max(pl.collect), "\n")
      }
      pl.tune = pl.collect[max.id]
    } else{ 
      lambda = lambda.former
      cat("----- Tuning END ----- \n")
      cat("Tuned lambda:", paste0(lambda, collapse = ", "), "\n")
      cat("Current pairwise log-likelihood:", pl.tune, "\n")
      cat("Maximum tried pairwise log-likelihood:", max(pl.collect), "\n")
      tuning = FALSE 
    }
  }
  
  if (!cv) {
    est.tune.trans = trans(est.tune, R.in, lambda)
  } else if (cv) {
    est.tune.trans = lapply(1:fold_n, function(i) trans(est.tune[[i]], R.in, lambda))
  }
  
  return(list(lambda = lambda, lambda.in = lambda.in,
              pl.tune = pl.tune, est.tune = est.tune, est.tune.trans = est.tune.trans,
              time = (proc.time() - ptm)[3]))
}
