Alg1 = function(R.in, reg.idx.in, mod.ini, mod.base, loc, sets, 
                pen.type, fold_n, stp.len = 0, model = "MSP"){
  
  # R.in: (integer): number of the subregion (in the base partition)
  # reg.idx.in (vector D): partition with subregion indicies
  # mod.ini (list): the model for initial values
  # mod.base (list): the model fitted on the base partition (under the same setting)
  # loc (matrix D*2): coordinate of the sites
  # sets (list): the setting of traning set, holdout set, validation set, and so on
  # pen.type (string): penalty used, "L1", "L2", or "None"
  # fold_n (integer): number of folds for cross-validation-type holdout sets
  # stp.len (numeric): between 0 and 1. The former corresponds to the self-adapted step length
  # model (string): "MSP" or "IMSP"
  
  ptm = proc.time()
  print("Merge subregions")
  
  est.ini = mod.ini$est.ini
  reg.idx.ini = mod.ini$reg.idx.ini
  obs.pairs = mod.ini$samplePairs.ns
  
  # previous
  R0 = R.in
  reg.idx0 = reg.idx.in
  adj.set0 = Neighbor(loc, reg.idx0)$adj.set
  lambda0 = mod.base$lambda
  lambda.former = mod.base$lambda.in
  
  
  # Update descending grid
  Lambda.list = UpdGrid(lambda0, lambda.former, Lambda.list)  
  
  # Current parameter estimates & penalized likelihood
  if (!cv) {
    est.ns0 = mod.base$est.tune
    est.ns0.trans = trans(est.ns0, R0, lambda0)
    ppl.ns0 = -PL.ns(est.ns0.trans, sets$D.hold, sets$Z.hold, reg.idx0[sets$holdIdx],
                     pen.type, adj.set0, lambda0, model = model)
    # pl.ns0 = mod.base$pl.tune
  } else if (cv) {
    est.ns0 = rowMeans(sapply(1:fold_n, function(i) mod.base$est.tune[[i]]))
    ppl.ns0 = mean(sapply(1:fold_n, function(i){
      -PL.ns(mod.base$est.tune.trans[[i]], sets$D.hold[[i]],
             sets$Z.hold[[i]], reg.idx0[sets$holdIdx[[i]]], pen.type, adj.set0,
             lambda0, model = model)
    }))
  }
  
  # Recording
  R.traj = c(R0)                        # Number of subregions
  reg.idx.traj = list(reg.idx0)         # Region indicator for each site
  adj.set.traj = list(adj.set0)         # Adjacent subregion pairs
  lambda.traj = list(lambda0)           # Optimal penalty parameter for accepted partitions
  est.ns.traj = list(est.ns0)           # Parameter estimates for accepted partitions
  ppl.ns.traj = c(ppl.ns0)              # PPL for accepted partitions
  
  
  # If KeepMerge = TRUE, we keep merging subregions, otherwise we stop
  KeepMerge = TRUE
  while (KeepMerge) {
    # Improvement = TRUE indicates that there is improvement in PPL, 
    # so we accept the candidate partition
    Improvement = FALSE
    threId = 1   # Indicator of threshold
    
    # Parameter estimate difference for neighboring subregions
    diff.pairs = sapply(1:nrow(adj.set0), function(i) {
      Diff(est.ns0, adj.set0[i,1], adj.set0[i,2], R0, lambda0) })
    
    diff.pairs.order = order(diff.pairs)
    diff.pairs.sort = sort(diff.pairs)
    
    # Merging thresholds
    if (stp.len == 0) {
      if (length(diff.pairs.sort) > 2) {
        clu <- kmeans(diff.pairs.sort, 2, nstart = 25, iter.max = 1000)
        sel.thre = which(clu$cluster == clu$cluster[1])
      } else {
        sel.thre = 1:length(diff.pairs.sort)
      }
      cat("\nSelf-adapted step len: ", 
          length(sel.thre)/length(diff.pairs.sort), '\n')
    } else {
      sel.thre = 1:(stp.len*length(diff.pairs.sort))
    }
    merge.thre = sort(unique(diff.pairs.sort[sel.thre]), decreasing = T)
    # ---
    
    cat("\nNew loop. Merge thresholds: ", merge.thre, '\n')
    
    while (!Improvement) {
      # ---
      tryThre = merge.thre[threId]
      underThre = which(diff.pairs <= tryThre)
      # Sort thresholds from large to smaller
      merge.idx = diff.pairs.order[1:length(underThre)]
      merge.pairs = adj.set0[merge.idx, ]
      
      # Merge subregions with given threshold
      merge.list = MergeSubr(est.ns0, merge.idx, merge.pairs, tryThre, reg.idx0, lambda0)
      reg.idx.try = merge.list$reg.idx
      R.try = length(unique(reg.idx.try))
      # ---
      cat("######## Try threshold", threId, ":", tryThre,
          "\n### Number of subregions: ", R.try, '\n')
      
      if (R.try == 1) {
        # Only one subregion left
        lambda.try = c(0, 0)
        est.ns.try = mod.ini$est.s
        est.ns.try.trans = trans(mod.ini$est.s)
        ppl.ns.try = mod.ini$pl.s
      } else {
        
        # If the partition has already been tested, pass
        # Some thresholds may cause exactly the same partition
        if (sum(reg.idx.try != reg.idx0) == 0) {
          cat("### Threshold", threId, ": the partition has been tested. \n")
          threId = threId + 1
          next }
        # ---
        
        # Fit nonstationary models for a candidate partition
        adj.set.try = Neighbor(loc, reg.idx.try)$adj.set
        
        if (!cv) {
          ini.val = Get_inival(est.ini, reg.idx.try, reg.idx.ini, lambda0)
          ini.val.trans = trans(ini.val, R.try, lambda0)
          
          mod.try = Fit.ns(ini.val, sets$D.train.ns, sets$Z.train.ns, 
                           R.try, reg.idx.try[sets$trainIdx.ns],
                           adj.set.try, pen.type, obs.pairs, 
                           lambda = lambda0, model = model) 
          pl.try = -PL.ns(mod.try$est.trans, sets$D.hold, sets$Z.hold, reg.idx.try[sets$holdIdx], 
                          pen.type = "None", adj.set.try, lambda0, model = model)
          
          lamTun.try = lamTuning(lambda0, Lambda.list, loc, reg.idx.in = reg.idx.try, 
                                 est.in = mod.try$est, pl.in = pl.try, pen.type = pen.type, 
                                 sets = sets, mod.ini = model.ini, cv = cv, show = T, 
                                 ncores = ncores, model = model)
          
          # Optimal lambda, parameter estimation, and corresponding PPL
          lambda.try = lamTun.try$lambda
          est.ns.try = lamTun.try$est.tune
          est.ns.try.trans = lamTun.try$est.tune.trans
          ppl.ns.try = -PL.ns(est.ns.try.trans, sets$D.hold, sets$Z.hold, reg.idx.try[sets$holdIdx], 
                              pen.type, adj.set.try, lambda.try, model = model)
          
        } else if (cv) {
          Mod.try = lapply(1:fold_n, function(i) {
            ini.val = Get_inival(est.ini, reg.idx.try, reg.idx.ini, lambda0)
            ini.val.trans = trans(ini.val, R.try, lambda0)
            mod.try = Fit.ns(ini.val, sets$D.train.ns[[i]], sets$Z.train.ns[[i]], 
                             R.try, reg.idx.try[sets$trainIdx.ns[[i]]], adj.set.try, 
                             pen.type, obs.pairs[[i]], lambda = lambda0, model = model) 
            pl.try = -PL.ns(mod.try$est.trans, sets$D.hold[[i]], sets$Z.hold[[i]], 
                            reg.idx.try[sets$holdIdx[[i]]], pen.type = "None", 
                            adj.set.try, lambda0, model = model)
            list(est = mod.try$est, pl = pl.try)
          })
            
          
          
          lamTun.try = lamTuning(lambda0, Lambda.list, loc, reg.idx.in = reg.idx.try, 
                                 est.in = lapply(1:fold_n, function(i) Mod.try[[i]]$est), 
                                 pl.in = mean(sapply(1:fold_n, function(i) Mod.try[[i]]$pl)), 
                                 pen.type = pen.type, sets = sets, mod.ini = model.ini, 
                                 cv = cv, show = T, ncores = ncores, model = model)
          
          # Optimal lambda, parameter estimation, and corresponding PPL
          lambda.try = lamTun.try$lambda
          est.ns.try = lamTun.try$est.tune
          est.ns.try.trans = lamTun.try$est.tune.trans
          ppl.ns.try = mean(sapply(1:fold_n, function(i) 
            -PL.ns(est.ns.try.trans[[i]], sets$D.hold[[i]], sets$Z.hold[[i]], 
                   reg.idx.try[sets$holdIdx[[i]]], pen.type, adj.set.try, 
                   lambda.try, model = model) ))
          
        }
        
        
      }

            
      # ------------------------------------------------------------------------
      # If there is improvement in PPL
      cat("###", threId, "/", length(merge.thre), '\n')
      print(threId == length(merge.thre))
      if (ppl.ns.try > ppl.ns0){ 
        cat('#############################################################', 
            "\nNew pl:", ppl.ns.try, "Old pl:", ppl.ns0, 
            '\nMerge subregions with threshold', merge.thre[threId], 
            '\nNumber of subregions:', R.try,
            '\n#############################################################')
        
        # Plot accepted partition
        plot(loc, col = sample(colors(), R.base)[reg.idx.try], pch=15)
        
        # Updation
        R0 = length(unique(reg.idx.try)) #
        reg.idx0 = sapply(reg.idx.try, function(id) which(unique(reg.idx0) == id))
        adj.set0 = Neighbor(loc, reg.idx0)$adj.set
        
        lambda.tmp = lambda0
        lambda0 = lambda.try
        Lambda.list = UpdGrid(lambda0, lambda.tmp, Lambda.list)
        
        if (!cv) {
          est.ns0 = est.ns.try
          ppl.ns0 = -PL.ns(est.ns.try.trans, sets$D.hold, sets$Z.hold, 
                           reg.idx0[sets$holdIdx], pen.type, adj.set0, lambda0, 
                           model = model)
        } else if (cv) {
          est.ns0 = rowMeans(sapply(1:fold_n, function(i) est.ns.try[[i]])) 
          ppl.ns0 = mean(sapply(1:fold_n, function(i){
            -PL.ns(est.ns.try.trans[[i]], sets$D.hold[[i]],
                   sets$Z.hold[[i]], reg.idx0[sets$holdIdx[[i]]], pen.type, adj.set0,
                   lambda0, model = model)
          }))
        }
        
        # Recording
        R.traj = c(R.traj, R0)
        reg.idx.traj = c(reg.idx.traj, list(reg.idx0))
        adj.set.traj = c(adj.set.traj, list(adj.set0))
        lambda.traj = c(lambda.traj, list(lambda0))
        est.ns.traj = c(est.ns.traj, list(est.ns0))
        ppl.ns.traj = c(ppl.ns.traj, ppl.ns0)
        
        Improvement = TRUE
        
        # Only one subregion left
        if (R0 == 1) {
          
        }
        
        
      } else if (threId == length(merge.thre)) { 
        # If no improvement for all thresholds, terminate.
        print('The merge procedure is terminated. A candidate is generated.\n')
        KeepMerge = FALSE
        break
      } else { 
        cat("Attempted ppl:", ppl.ns.try, "Current ppl:", ppl.ns0, '\n')
        cat("Reject the merge with threshold", merge.thre[threId], '\n')
        cat("=> Threshold", threId, "failed. \n")
        threId = threId + 1 
        # cat("###", threId, "/", length(merge.thre), '\n')
      }
      
    }
  }
  
  # stopCluster(cl)
  RT = (proc.time() - ptm)[3]
  print(RT)
  
  list(R.traj = R.traj, 
       reg.idx.traj = reg.idx.traj,
       adj.set.traj = adj.set.traj,
       lambda.traj = lambda.traj,
       est.ns.traj = est.ns.traj,
       ppl.ns.traj = ppl.ns.traj,
       stp.len = stp.len,
       RT = RT)
}



