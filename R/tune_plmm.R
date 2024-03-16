tune_plmm <- function(data, gamma, lambda, crit = "BIC", intercept = T, timexgroup = T) {
  grid.param <- expand.grid(lambda, gamma)
  list.MixteNonPara <- list()
  
  Dico.norm <- CreationBases(data$position)
  F.Bases.norm <- Dico.norm$F.Bases
  
  for (i in 1:nrow(grid.param)) {
    lambda.i <- grid.param[i, ][[1]]
    gamma.i <- grid.param[i, ][[2]]
    
    # start_time <- Sys.time()
    
    list.MixteNonPara[[i]] <- plmmlasso(
      Y = data$Y, series = data$series, position = data$position,
      X = subset(data, select = -c(Y, series, position)),
      F.Bases = F.Bases.norm, gam.cste = gamma.i, intercept = intercept,
      lambda.grid = lambda.i, timexgroup = timexgroup
    )
    
    # end_time <- Sys.time()
    
    # list.MixteNonPara[[i]]$time <- difftime(end_time, start_time, units = "secs")
  }
  
  list.MixteNonPara <- list.MixteNonPara[unlist(lapply(
    list.MixteNonPara,
    function(x) x$converged
  ))]
  
  best.lpmm <- which.min(unlist(lapply(list.MixteNonPara, function(x) x[crit])))
  
  return(list.MixteNonPara[[best.lpmm]])
}
