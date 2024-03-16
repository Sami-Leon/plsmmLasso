calc_criterion <- function(name, out.EM, data, logLik, nonpara = F) {
  ntot <- nrow(data)
  k <- sum(out.EM$theta != 0)
  d <- length(out.EM$theta)
  
  if (nonpara == T) {
    k <- k + length(out.EM$FunctSelect)
    d <- length(out.EM$theta) + length(as.vector(out.EM$Coef.Val))
  }
  
  if (name == "BIC") {
    return(-2 * logLik + log(ntot) * k)
  }
  if (name == "BICC") {
    return(-2 * logLik + max(1, log(log(d))) * log(ntot) * k)
  }
  if (name == "EBIC") {
    return(-2 * logLik + (log(ntot) + 2 * log(d)) * k)
  }
}

init_params <- function(Data) {
  Yi.dot <- tapply(Data$Y,
                   list(series = as.factor(Data$series)),
                   FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  su <- var(c(Yi.dot), na.rm = T)
  se <- var(Data$Y, na.rm = T) - su
  g <- se / su
  invisible(list(g = g, se = se))
}

E_step <- function(Data, out.F, g, ni, theta, fixed_effects) {
  X <- as.matrix(cbind(1, subset(Data, select = -c(series, position, Y))))
  
  if (fixed_effects == TRUE) {
    res <- data.frame(
      series = Data$series,
      position = Data$position,
      resid = (Data$Y - out.F$F.fit - X %*% theta)
    )
  } else {
    res <- data.frame(
      series = Data$series,
      position = Data$position,
      resid = (Data$Y - out.F$F.fit)
    )
  }
  
  
  U2 <- tapply(res$resid,
               list(series = as.factor(res$series)),
               FUN = function(x) sum(x, na.rm = TRUE)
  ) / (ni + g)
  
  U2 <- data.frame(
    series = unique(Data$series),
    U2 = c(t(U2))
  )
  
  invisible(list(U2 = U2))
}

M_step_standard_error <- function(Data, out.F, g, se, U2, ni, theta, fixed_effects) {
  X <- as.matrix(cbind(1, subset(Data, select = -c(series, position, Y))))
  
  I <- length(unique(Data$series))
  N <- sum(ni)
  
  if (fixed_effects == TRUE) {
    res <- data.frame(
      series = Data$series,
      position = Data$position,
      resid = (Data$Y - out.F$F.fit - X %*% theta)
    )
  } else {
    res <- data.frame(
      series = Data$series,
      position = Data$position,
      resid = (Data$Y - out.F$F.fit)
    )
  }
  
  repU2 <- rep(U2$U2, ni)
  
  se <- sum((res$resid - repU2)^2, na.rm = TRUE) + sum((ni * se) / (ni + g))
  se <- se / N
  
  se[is.na(se)] <- 0
  se[se == Inf] <- 0
  
  return(se = se)
}

M_step_random_effects <- function(Data, g, se, U2, ni) {
  I <- length(unique(Data$series))
  su <- sum((U2$U2)^2, na.rm = T) + sum(se / (ni + g), na.rm = TRUE)
  su <- su / I
  
  su[su == Inf] <- 0
  return(su = su)
}

offset_random_effects <- function(Data, U2, ni, fixed_effects) {
  X <- as.matrix(cbind(1, subset(Data, select = -c(series, position, Y))))
  Data.tmp <- Data
  U2.tmp <- rep(U2$U2, ni)
  
  if (fixed_effects == TRUE) {
    Data.tmp$Y <- Data.tmp$Y - U2.tmp
  } else {
    Data.tmp$Y <- Data.tmp$Y - U2.tmp
  }
  
  return(Data.tmp = Data.tmp)
}

joint_lasso <- function(Data, F.Bases, se, gam.cste, intercept, fixed_effects,
                      lambda.grid, timexgroup, pre.D, debias = FALSE) {
  M <- (dim(F.Bases)[2]) / 2
  
  X <- as.matrix(subset(Data, select = -c(series, position, Y)))
  
  X.stand <- scale(X, scale = T)
  
  X.mean <- attr(X.stand, "scaled:center")
  X.sd <- attr(X.stand, "scaled:scale")
  
  X <- cbind(1, X)
  X.stand <- cbind(1, X.stand)
  
  Flars.X <- cbind(X.stand, F.Bases)
  
  D <- diag(1, nrow = ncol(Flars.X), ncol = ncol(Flars.X))
  
  Lambda <- sqrt(se * gam.cste * log(M))
  
  D[(ncol(X.stand) + 1):ncol(D), (ncol(X.stand) + 1):ncol(D)] <- pre.D * (Lambda / lambda.grid)
  
  D.inv <- D
  diag(D.inv) <- 1 / diag(D.inv)
  
  Flars.X.lasso <- Flars.X %*% D.inv
  
  y.stand <- scale(Data$Y)
  y.mean <- attr(y.stand, "scaled:center")
  y.sd <- attr(y.stand, "scaled:scale")
  
  if (intercept == TRUE) {
    coef.glmnet.w <- as.vector(coef(glmnet::glmnet(Flars.X.lasso[, -1],
                                                   y.stand,
                                                   alpha = 1, lambda = lambda.grid,
                                                   standardize = FALSE,
                                                   intercept = TRUE
    )))
  } else {
    coef.glmnet.w <- as.vector(coef(glmnet::glmnet(Flars.X.lasso[, -1],
                                                   y.stand,
                                                   alpha = 1, lambda = lambda.grid,
                                                   standardize = FALSE,
                                                   intercept = FALSE
    )))
  }
  
  coef.glmnet.w[1] <- (coef.glmnet.w[1] - sum((X.mean / X.sd) * coef.glmnet.w[2:ncol(X)])) * y.sd + y.mean
  
  coef.lasso <- D.inv %*% coef.glmnet.w
  coef.lasso[-1] <- coef.lasso[-1] * y.sd
  coef.lasso[2:ncol(X), 1] <- coef.lasso[2:ncol(X), 1] / X.sd
  
  theta <- coef.lasso[1:ncol(X), 1]
  names(theta) <- colnames(X)
  names(theta)[1] <- "Intercept"
  
  Coef.Val <- coef.lasso[(ncol(X) + 1):nrow(coef.lasso), 1]
  
  theta[-1][abs(theta)[-1] < 0.01] <- 0
  
  F.esti.lars.w <- F.Bases %*% Coef.Val
  
  out.F <- data.frame(
    position = Data$position,
    F.fit = F.esti.lars.w,
    Group = Data[, "Group"]
  )
  
  X.fit <- X %*% theta
  
  FunctSelect <- which(Coef.Val != 0)
  
  invisible(list(
    out.F = out.F, FunctSelect = FunctSelect, Coef.Val = Coef.Val,
    theta = theta, X.fit = X.fit
  ))
}

plmm_lasso <- function(Y, series, position, X = NULL, F.Bases, gam.cste, intercept,
                      lambda.grid, timexgroup, tol.EM = 0.001) {
  fixed_effects <- ifelse(is.null(X), FALSE, TRUE)
  
  Data <- as.data.frame(cbind(Y, series, position, X))
  ni <- as.vector(table(Data$series))
  
  if (timexgroup) {
    F.Bases.timexgroup <- matrix(nrow = nrow(Data), ncol = ncol(F.Bases) * 2)
    for (i in 1:nrow(Data)) {
      if (Data$Group[i] == Data$Group[1]) {
        F.Bases.timexgroup[i, ] <- c(F.Bases[i, ], rep(0, ncol(F.Bases)))
      } else {
        F.Bases.timexgroup[i, ] <- c(rep(0, ncol(F.Bases)), F.Bases[i, ])
      }
    }
  } else {
    F.Bases.timexgroup <- F.Bases
  }
  
  pre.D <- diag(sqrt(apply(F.Bases.timexgroup^2, 2, sum)))
  
  ## Initialization
  out.si <- init_params(Data)
  g <- out.si[[1]]
  se <- out.si[[2]]
  
  theta <- rep(0, dim(Data)[2] - 2)
  
  init.EstiF <- joint_lasso(
    Data = Data, F.Bases = F.Bases.timexgroup, se = se,
    gam.cste = gam.cste, intercept = intercept,
    fixed_effects = fixed_effects, lambda.grid = lambda.grid,
    timexgroup = timexgroup, pre.D = pre.D
  )
  
  out.F <- init.EstiF$out.F
  theta <- init.EstiF$theta
  
  
  maxIter <- 50
  delta.EM <- Inf
  Iter <- 0
  
  theta.MSE <- NULL
  F.fit.MSE <- NULL
  U2.MSE <- NULL
  overall.MSE <- NULL
  
  delta.EM.iter <- NULL
  
  while ((delta.EM > tol.EM) & (Iter < maxIter)) {
    Iter <- Iter + 1
    
    out.E <- E_step(Data, out.F, g, ni, theta, fixed_effects)
    U2 <- out.E$U2
    
    su.tmp <- M_step_random_effects(Data, g, se, U2, ni)
    
    se.tmp <- M_step_standard_error(
      Data, out.F, g, se, U2, ni, theta,
      fixed_effects
    )
    g.tmp <- se.tmp / su.tmp
    g.tmp[is.nan(g.tmp)] <- 0
    su.tmp[is.nan(su.tmp)] <- 0
    
    # F estimation
    out.Data <- offset_random_effects(Data, U2, ni, fixed_effects)
    
    Res.F <- joint_lasso(
      Data = out.Data, F.Bases = F.Bases.timexgroup, se = se.tmp, gam.cste = gam.cste,
      intercept = intercept, fixed_effects = fixed_effects,
      lambda.grid = lambda.grid, timexgroup = timexgroup,
      pre.D = pre.D
    )
    
    
    out.F.tmp <- Res.F$out.F
    
    theta.tmp <- Res.F$theta
    
    
    delta.F <- 0
    delta.theta <- 0
    delta.se <- 0
    delta.su <- 0
    
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
      delta.EM <- sum((tp0 - tp1)^2)
      
      delta.F <- sum((out.F$F.fit - out.F.tmp$F.fit)^2)
      delta.theta <- sum((theta - theta.tmp)^2)
      delta.se <- sum((se - se.tmp)^2)
      delta.su <- sum((se / g - se.tmp / g.tmp)^2)
    }
    
    if ((delta.EM > 3) & (Iter > 3)) {
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
    
    cat("Iter ", Iter, delta.EM, delta.F, delta.theta, delta.se, delta.su, "\n")
    
    delta.EM.iter[Iter] <- delta.EM
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
  converged <- ifelse(Iter >= maxIter, F, T)
  
  Z <- model.matrix(~ 0 + factor(series), Data)
  logLik <- mvtnorm::dmvnorm(
    x = Data$Y,
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
    Res.F = Res.F, se = se, su = su, U2 = U2, ni = ni, delta.EM = delta.EM.iter,
    hyper.parameters = hyper.parameters, converged = converged, BIC = BIC,
    BIC.nonpara = BIC.nonpara, BICC = BICC, BICC.nonpara = BICC.nonpara,
    EBIC = EBIC, EBIC.nonpara = EBIC.nonpara
  ))
}
