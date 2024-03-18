calc_criterion <- function(crit, res_output, log_lik, nonpara = FALSE) {
  n <- res_output$X.fit
  k <- sum(res_output$theta != 0)
  d <- length(res_output$theta)
  
  if (nonpara == TRUE) {
    k <- k + length(res_output$FunctSelect)
    d <- length(res_output$theta) + length(as.vector(res_output$Coef.Val))
  }
  
  if (crit == "BIC") {
    return(-2 * log_lik + log(n) * k)
  }
  if (crit == "BICC") {
    return(-2 * log_lik + max(1, log(log(d))) * log(n) * k)
  }
  if (crit == "EBIC") {
    return(-2 * log_lik + (log(n) + 2 * log(d)) * k)
  }
}

init_params <- function(y, series) {
  y_series_means <- tapply(y,
                   list(series = as.factor(series)),
                   FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  su <- var(c(y_series_means), na.rm = T)
  se <- var(y, na.rm = T) - su
  sr <- se / su
  invisible(list(g = sr, se = se))
}

E_step <- function(x, y, series, out.F, g, ni, theta) {
  
  res <- data.frame(
      series = series,
      resid = (y - out.F$F.fit - x %*% theta)
    )

  U2 <- tapply(res$resid,
               list(series = as.factor(res$series)),
               FUN = function(x) sum(x, na.rm = TRUE)
  ) / (ni + g)
  
  U2 <- data.frame(
    series = unique(Data$series),
    U2 = c(t(U2))
  )
  
  invisible(U2)
}

M_step_standard_error <- function(x, y, out.F, g, se, U2, ni, theta) {
  n <- length(y)
  
  repU2 <- rep(U2$U2, ni)
  
  se <- sum((y - out.F$F.fit - x %*% theta - repU2)^2, na.rm = TRUE) + sum((ni * se) / (ni + g))
  se <- se / n
  
  se[is.na(se) | se == Inf] <- 0
  
  return(se = se)
}

M_step_random_effects <- function(series, g, se, U2, ni) {
  N <- length(unique(series))
  su <- sum((U2$U2)^2, na.rm = T) + sum(se / (ni + g), na.rm = TRUE)
  su <- su / N
  
  su[is.na(se) | se == Inf] <- 0
  return(su = su)
}

offset_random_effects <- function(y, U2, ni) {
  repU2 <- rep(U2$U2, ni)
  y <- y - repU2
  return(y)
}

joint_lasso <- function(x, y, t, name_group_var, bases, se, gamma,
                      lambda, pre_D) {

  x_stand <- scale(x, scale = TRUE)
  
  x_mean <- attr(x_stand, "scaled:center")
  x_sd <- attr(x_stand, "scaled:scale")
  
  x <- cbind(1, x)
  x_stand <- cbind(1, x_stand)

  combined_x_bases <- cbind(x_stand, bases)
  
  p = ncol(x_stand)
  M <- ncol(bases)
  pM = p + M 
  
  D <- diag(1, nrow = pM, ncol = pM)
  
  D[(p + 1):pM, (p + 1):pM] <- pre_D * (sqrt(se * gamma * log(M/2)) / lambda)
  
  D_inv <- D
  diag(D_inv) <- 1 / diag(D_inv)
  
  combined_x_bases_lasso <- combined_x_bases %*% D_inv
  
  y_stand <- scale(y)
  y_mean <- attr(y_stand, "scaled:center")
  y_sd <- attr(y_stand, "scaled:scale")
  
  coef_joint_lasso <- as.vector(coef(glmnet::glmnet(combined_x_bases_lasso[, -1],
                                                   y_stand,
                                                   alpha = 1, lambda = lambda,
                                                   standardize = FALSE,
                                                   intercept = TRUE
    )))

  
  coef_joint_lasso[1] <- (coef_joint_lasso[1] - sum((x_mean / x_sd) * coef_joint_lasso[2:p])) * y_sd + y_mean
  
  coef_joint_lasso <- D_inv %*% coef_joint_lasso
  coef_joint_lasso[-1] <- coef_joint_lasso[-1] * y_sd
  coef_joint_lasso[2:p, 1] <- coef_joint_lasso[2:p, 1] / x_sd
  
  theta <- coef_joint_lasso[1:p, 1]
  names(theta) <- colnames(x)
  names(theta)[1] <- "Intercept"
  
  alpha <- coef_joint_lasso[(p + 1):nrow(coef_joint_lasso), 1]
  
  # theta[-1][abs(theta)[-1] < 0.01] <- 0
  
  f_hat <- bases %*% alpha
  
  out.F <- data.frame(
    position = t,
    F.fit = f_hat,
    Group = x[, name_group_var]
  )
  
  x_fit <- x %*% theta
  
  selected_functions <- which(alpha != 0)
  
  invisible(list(
    out.F = out.F, FunctSelect = selected_functions, Coef.Val = alpha,
    theta = theta, X.fit = x_fit
  ))
}

plmm_lasso <- function(x, y, series, t, name_group_var, bases, 
                       gamma, lambda, timexgroup, tol.EM = 0.001) {
  
  # Data <- as.data.frame(cbind(y, series, position, X))
  p = ncol(x)
  ni <- as.vector(table(series))

  if (timexgroup) {
    n = length(y)
    vec_group =  x[,name_group_var]
    ref_group = vec_group[1]
    ncol(bases) = M
    bases_timexgroup <- matrix(nrow = n, ncol = M * 2)
    for (i in 1:n) {
      if (vec_group[i] == ref_group) {
        bases_timexgroup[i, ] <- c(bases[i, ], rep(0, M))
      } else {
        bases_timexgroup[i, ] <- c(rep(0, M), bases[i, ])
      }
    }
  }
  
  bases = bases_timexgroup

  
  pre_D <- diag(sqrt(apply(bases^2, 2, sum)))
  
  ## Initialization
  out.si <- init_params(y = y, series = series)
  g <- out.si[[1]]
  se <- out.si[[2]]
  
  theta <- rep(0, p + 1)
  
  init.EstiF <- joint_lasso(x = x, y = y, t = t, name_group_var = name_group_var, 
                            bases = bases, se = se, gamma = gamma,lambda = lambda,
                            pre_D = pre_D)
  
  out.F <- init.EstiF$out.F
  theta <- init.EstiF$theta
  
  
  max_iter <- 50
  cvg_crit <- Inf
  Iter <- 0
  
  while ((cvg_crit > tol.EM) & (Iter < max_iter)) {
    Iter <- Iter + 1
    
    out.E <- E_step(x = x, y = y, series = series, out.F = out.F, g = g, ni = ni,
                    theta = theta)
    U2 <- out.E
    # here
    se.tmp <- M_step_standard_error(x = x, y = y, out.F = out.F, g = g, se = se,
                                    U2 = U2, ni = ni, theta = theta)
    
    su.tmp <- M_step_random_effects(series = series, g = g, se = se, U2 = U2, 
                                    ni = ni)
    
    g.tmp <- se.tmp / su.tmp
    # g.tmp[is.nan(g.tmp)] <- 0
    # su.tmp[is.nan(su.tmp)] <- 0
    
    # F estimation
    y_offset <- offset_random_effects(y = y, U2 = U2, ni = ni)
    
    Res.F <- joint_lasso(x = x, y = y_offset, t = t, name_group_var = name_group_var,
                         bases = bases, se = se, gamma = gamma, lambda = lambda, 
                         pre_D = pre_D)
    
    
    out.F.tmp <- Res.F$out.F
    
    theta.tmp <- Res.F$theta
    
    delta_f <- 0
    delta_theta <- 0
    delta_se <- 0
    delta_su <- 0
    
    if (Iter == 2) {
      t2 <- c(out.F$F.fit, se, se / g, theta)
    }
    if (Iter == 3) {
      t1 <- c(out.F$F.fit, se, se / g, theta)
      t0 <- c(out.F.tmp$F.fit, se.tmp, se.tmp / g.tmp, theta.tmp)
      tp0 <- (t2 - t1) / sum(((t2 - t1)^2)) + (t0 - t1) / sum(((t0 - t1)^2))
      tp0 <- t1 + tp0 / sum(tp0^2)
    }
    if (Iter > 3) {
      t2 <- t1
      t1 <- t0
      t0 <- c(out.F.tmp$F.fit, se.tmp, se.tmp / g.tmp, theta.tmp)
      tp1 <- tp0
      tp0 <- (t2 - t1) / sum(((t2 - t1)^2)) + (t0 - t1) / sum(((t0 - t1)^2))
      tp0 <- t1 + tp0 / sum(tp0^2)
      cvg_crit <- sum((tp0 - tp1)^2)
      
      delta_f <- sum((out.F$F.fit - out.F.tmp$F.fit)^2)
      delta_theta <- sum((theta - theta.tmp)^2)
      delta_se <- sum((se - se.tmp)^2)
      delta_su <- sum((se / g - se.tmp / g.tmp)^2)
    }
    
    if ((cvg_crit > 3) & (Iter > 3)) {
      Iter <- Iter + (50 - Iter)
    }
    
    # Update
    se <- se.tmp
    su <- su.tmp
    out.F <- out.F.tmp
    g <- g.tmp
    theta <- theta.tmp
    
    mean.0.tmp <- attr(scale(unique(out.F[out.F$Group == 0, ]$F.fit),
                             scale = F
    ), "scaled:center")
    mean.1.tmp <- attr(scale(unique(out.F[out.F$Group == 1, ]$F.fit),
                             scale = F
    ), "scaled:center")
    
    out.F.tmp <- out.F
    
    out.F.tmp[out.F.tmp$Group == 0, ]$F.fit <- out.F.tmp[out.F.tmp$Group == 0, ]$F.fit - mean.0.tmp
    out.F.tmp[out.F.tmp$Group == 1, ]$F.fit <- out.F.tmp[out.F.tmp$Group == 1, ]$F.fit - mean.1.tmp
    
    theta.tmp <- theta
    
    theta.tmp[2] <- theta.tmp[2] + (mean.1.tmp - mean.0.tmp)
    theta.tmp[1] <- theta.tmp[1] + mean.0.tmp
    
    cat("Iter ", Iter, cvg_crit, delta_f, delta_theta, delta_se, delta_su, "\n")
  }
  
  mean.0 <- attr(scale(unique(Res.F$out.F[Res.F$out.F$Group == 0, ]$F.fit),
                       scale = F
  ), "scaled:center")
  mean.1 <- attr(scale(unique(Res.F$out.F[Res.F$out.F$Group == 1, ]$F.fit),
                       scale = F
  ), "scaled:center")
  
  Res.F$out.F[Res.F$out.F$Group == 0, ]$F.fit <- Res.F$out.F[Res.F$out.F$Group == 0, ]$F.fit - mean.0
  Res.F$out.F[Res.F$out.F$Group == 1, ]$F.fit <- Res.F$out.F[Res.F$out.F$Group == 1, ]$F.fit - mean.1
  
  Res.F$theta["Group"] <- Res.F$theta["Group"] + (mean.1 - mean.0)
  Res.F$theta["Intercept"] <- Res.F$theta["Intercept"] + mean.0
  Res.F$X.fit <- as.matrix(cbind(1, X)) %*% Res.F$theta
  
  hyper.parameters <- data.frame(lambda.grid = lambda.grid, gam.cste = gam.cste)
  converged <- ifelse(Iter >= max_iter, F, T)
  
  Z <- model.matrix(~ 0 + factor(series), Data)
  logLik <- mvtnorm::dmvnorm(
    x = y,
    mean = as.vector(Res.F$X.fit) + Res.F$out.F$F.fit,
    sigma = diag(nrow(Z)) * se + su * Z %*% t(Z), log = T
  )
  
  BIC <- calc_criterion(
    name = "BIC", out.EM = Res.F, data = Data,
    logLik = logLik, nonpara = F
  )
  BIC.nonpara <- calc_criterion(
    name = "BIC", out.EM = Res.F, data = Data,
    logLik = logLik, nonpara = T
  )
  BICC <- calc_criterion(
    name = "BICC", out.EM = Res.F, data = Data,
    logLik = logLik, nonpara = F
  )
  BICC.nonpara <- calc_criterion(
    name = "BICC", out.EM = Res.F, data = Data,
    logLik = logLik, nonpara = T
  )
  EBIC <- calc_criterion(
    name = "EBIC", out.EM = Res.F, data = Data, logLik =
      logLik, nonpara = F
  )
  EBIC.nonpara <- calc_criterion(
    name = "EBIC", out.EM = Res.F, data = Data,
    logLik = logLik, nonpara = T
  )
  
  return(list(
    Res.F = Res.F, se = se, su = su, U2 = U2, ni = ni,
    hyper.parameters = hyper.parameters, converged = converged, BIC = BIC,
    BIC.nonpara = BIC.nonpara, BICC = BICC, BICC.nonpara = BICC.nonpara,
    EBIC = EBIC, EBIC.nonpara = EBIC.nonpara
  ))
}
