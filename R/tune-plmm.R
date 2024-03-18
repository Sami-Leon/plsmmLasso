tune_plmm <- function(x, y, series, t, name_group_var, bases,
                      gamma_vec, lambda_vec, timexgroup, criterion, ...) {
  grid_hyper <- expand.grid(lambda_vec, gamma_vec)
  plmm_models <- vector("list", nrow(grid_hyper))

  for (i in 1:nrow(grid_hyper)) {
    lambda_i <- grid_hyper[i, ][[1]]
    gamma_i <- grid_hyper[i, ][[2]]


    plmm_models[[i]] <- plmm_lasso(
      x = x, y = y, series = series, t = t,
      name_group_var = name_group_var, bases = bases,
      gamma = gamma_i, lambda = lambda_i,
      timexgroup = timexgroup, criterion = criterion, ...
    )
  }

  plmm_models <- plmm_models[unlist(lapply(
    plmm_models,
    function(x) x$converged
  ))]

  best_plmm <- which.min(unlist(lapply(plmm_models, function(x) x[criterion])))

  return(plmm_models[[best_plmm]])
}
