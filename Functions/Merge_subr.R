UpdGrid = function(lambda, lambda.former, Lambda.list){
  # lambda (vector): optimal lambda after lambda tuning
  # lambda.former (vector): lambda before lambda tuning
  # Lambda.list (list): descending grid
  
  # If lambda_i is different from lambda.former_i, descending grid should be updated
  for (i in length(lambda)){
    i = 1
    if (lambda[i] != lambda.former[i]){
      grid.pos = which(Lambda.list[[i]] == lambda[i])
      if (grid.pos != length(Lambda.list[[i]])) { #
        lower.grid = Lambda.list[[i]][grid.pos+1]
        Lambda.list[[i]] = c(Lambda.list[[i]][1:grid.pos], (lambda[i]+lower.grid)/2, 
                             Lambda.list[[i]][(grid.pos+1):length(Lambda.list[[i]])])
      }
    }
  }
  return(Lambda.list)
}

diff.scale = function(vec, diff.type = "log"){ 
  if (diff.type == "identity") { vec
  } else if (diff.type == "log") { log(vec)
  }
}

# Measure estimate difference of neighboring subregions
Diff = function(par, r1, r2, R.in, lambda, diff.type = "log"){
  # par: parameter estimates two subregions
  # r1 & r2 (number): subregion 1 and subregion 2
  # R.in (number): number of subregions
  # diff.type (string): dufference type, "identity" or "log". Others are possible.
  
  if (lambda[1] == Inf){
    diff1 = 0; pos = 1
  } else if (lambda[1] != Inf) {
    diff1 = abs(diff.scale(par[r1], diff.type) -
                  diff.scale(par[r2], diff.type)); pos = R.in}
  if (lambda[2] == Inf){
    diff2 = 0
  } else if (lambda[2] != Inf) {
    diff2 = abs(diff.scale(par[pos + r1], diff.type) -
                  diff.scale(par[pos + r2], diff.type)) }
  diff1 + diff2
}

MergeSubr = function(est.in, merge.idx, merge.pairs, thre, reg.idx.in, lambda, diff.type = "log"){
  
  # est.in (vector): parameter estimates for each subregion
  # merge.idx (vector): indices of subregion pairs that will be merged
  # merge.pairs (matrix with 2 columns): subregion pairs that will be merged
  # thre (number): merging threshold
  # reg.idx.in (vector D): current partition
  
  reg.idx.in.uni = unique(reg.idx.in)
  R.in = length(reg.idx.in.uni)
  
  # ---
  est.mat = vec2mat(est.in, R.in, lambda)
  
  grouped = c(); groups = list()
  
  for (r in 1:length(merge.idx)){ 
    merge.pair = as.matrix(merge.pairs[r, ])
    
    r1 = merge.pair[1]; r2 = merge.pair[2] #pos
    if (!(r1 %in% grouped) & !(r2 %in% grouped)) {
      est.tmp = est.mat[,merge.pair]
      groups = c(groups, list(list(gp = merge.pair, est.tmp = est.tmp)))
      grouped = c(grouped, merge.pair)
    } else if (r1 %in% grouped & !(r2 %in% grouped)) {
      gId = which(sapply(1:length(groups), function(i) r1 %in% groups[[i]]$gp))
      d2Group = sapply(groups[[gId]]$gp, function(ri) { Diff(est.in, r2, ri, R.in, lambda, diff.type) })
      if (max(d2Group) <= thre) {
        groups[[gId]]$gp = c(groups[[gId]]$gp, r2)
        groups[[gId]]$est.tmp = cbind(groups[[gId]]$est.tmp, est.mat[,r2])
        grouped = c(grouped, r2)
      }
    } else if (!(r1 %in% grouped) & r2 %in% grouped) {
      gId = which(sapply(1:length(groups), function(i) r2 %in% groups[[i]]$gp))
      d2Group = sapply(groups[[gId]]$gp, function(ri) { Diff(est.in, r1, ri, R.in, lambda, diff.type) })
      if (max(d2Group) <= thre) {
        groups[[gId]]$gp = c(groups[[gId]]$gp, r1)
        groups[[gId]]$est.tmp = cbind(groups[[gId]]$est.tmp, est.mat[,r1])
        grouped = c(grouped, r1)
      }
    } else if (r1 %in% grouped & r2 %in% grouped) {
      gId1 = which(sapply(1:length(groups), function(i) r1 %in% groups[[i]]$gp))
      gId2 = which(sapply(1:length(groups), function(i) r2 %in% groups[[i]]$gp))
      if (gId1 != gId2) {
        d2Group = sapply(groups[[gId1]]$gp, function(ri) {
          sapply(groups[[gId2]]$gp, function(rj) { Diff(est.in, ri, rj, R.in, lambda, diff.type) }) })
        if (max(d2Group) <= thre) {
          groups[[gId1]]$gp = c(groups[[gId1]]$gp, groups[[gId2]]$gp)
          groups[[gId1]]$est.tmp = cbind(groups[[gId1]]$est.tmp, groups[[gId2]]$est.tmp)
          groups[[gId2]] = NULL
        }
      }
    }
  }
  
  if (length(grouped) != R.in){
    for (r in reg.idx.in.uni[-grouped]){
      est.tmp = est.mat[,r]
      groups = c(groups, list(list(gp = r, est.tmp = est.tmp)))
      grouped = c(grouped, r)
    }
  }
  
  reg.idx.tmp = rep(NaN, length(reg.idx.in))
  for (r in 1:length(groups)){
    reg.idx.tmp[which(reg.idx.in %in% groups[[r]]$gp)] = r
  }
  
  reg.idx.tmp = sapply(reg.idx.tmp, function(id) which(unique(reg.idx.tmp) == id))
  
  list(reg.idx = reg.idx.tmp, grouped = grouped, groups = groups)
}



