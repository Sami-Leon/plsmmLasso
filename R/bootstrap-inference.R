boot_samples <- function(data, n_boot) {
  unique_patients <- unique(data$series)
  n_series <- length(unique_series)
  
  bootstrap_samples <- lapply(1:n_boot, function(x) {
    sample_id <- sample(unique_patients, n_series, replace = TRUE)
    mat_boot <- do.call(rbind, lapply(sample_id, function(j) data[data$series == j, ]))
    return(mat_boot)
  })
  
  return(bootstrap_samples)
}

boot_pre <- function(data, timexgroup = TRUE) {
  bases <- create_bases(data$t)$bases
  
  if (timexgroup) {
    n <- nrow(data)
    vec_group <- data$group
    ref_group <- vec_group[1]
    M <- ncol(bases)
    
    # Create a logical index indicating where vec_group is equal to ref_group
    index_ref_group <- vec_group == ref_group
    
    # Initialize bases_timexgroup matrix
    bases_timexgroup <- matrix(0, nrow = n, ncol = M * 2)
    
    # Fill the bases_timexgroup matrix using vectorized operations
    bases_timexgroup[index_ref_group, 1:M] <- bases[index_ref_group, ]
    bases_timexgroup[!index_ref_group, (M+1):(2*M)] <- bases[!index_ref_group, ]
    
    bases <- bases_timexgroup
  }
  
  
  return(bases)
}


fit.boot <- function(Data, sig.init) {
  pre.boot <- boot_pre(Data)
  # sig.init<-scalreg::scalreg(scale(pre.boot), scale(Data$Y))$hsigma
  
  lambda.grid = seq(1.2, sig.init, -0.1)*sqrt(2*log(ncol(pre.boot))/length(Data$Y))
  
  cv_fit <- glmnet::cv.glmnet(pre.boot, Data$Y, alpha = 1, 
                              lambda = lambda.grid)
  
  final_fit <- glmnet::glmnet(pre.boot, Data$Y, alpha = 1, lambda = cv_fit$lambda.min)
  fitted_values <- glmnet::predict.glmnet(final_fit, newx = pre.boot, s = cv_fit$lambda.min)
  
  Data$F.fit <- fitted_values
  
  out.F <- Data[, c("position", "F.fit", "Group")]
  
  out.F[out.F$Group == 0, ]$F.fit <- out.F[out.F$Group == 0, ]$F.fit - attr(scale(unique(out.F[out.F$Group == 0, ]$F.fit),
                                                                                  scale = F
  ), "scaled:center")
  out.F[out.F$Group == 1, ]$F.fit <- out.F[out.F$Group == 1, ]$F.fit - attr(scale(unique(out.F[out.F$Group == 1, ]$F.fit),
                                                                                  scale = F
  ), "scaled:center")
  return(list(out.F = out.F, Coef.Val = as.vector(final_fit$beta[, 1])))
}

diff.f <- function(res.boot) {
  fit0 <- res.boot[res.boot$Group == 0, ]
  fit0 <- fit0[!duplicated(fit0$position), ]
  fit0 <- fit0[order(fit0$position), ]
  
  fit1 <- res.boot[res.boot$Group == 1, ]
  fit1 <- fit1[!duplicated(fit1$position), ]
  fit1 <- fit1[order(fit1$position), ]
  
  out <- cbind(fit1$position, diff = fit0$F.fit - fit1$F.fit)
  
  colnames(out) <- c("position", "diff")
  out <- as.data.frame(out)
  return(out)
}

overall.test <- function(list.boot.fit, model) {
  df_list_fit <- lapply(seq_along(list.boot.fit), function(i) {
    df <- list.boot.fit[[i]]$out.F
    df$boot <- as.character(i) # Create a group variable
    return(df)
  })
  
  diff.list <- lapply(df_list_fit, diff.f)
  
  Tobs <- sum((diff.f(model$Res.F$out.F)$diff)^2)
  Tboot <- unlist(lapply(diff.list, function(x) {
    sum(x$diff^2)
  }))
  
  pvalue <- pnorm(abs((2 * Tobs - mean(Tboot)) / sqrt(var(Tboot))),
                  lower.tail = F
  ) * 2
  
  df.overall.f <- data.frame(Tobs, pvalue)
  
  colnames(df.overall.f) <- c("T", "p-value")
  
  return(df.overall.f)
}

pred.f <- function(model, data, byseq = 0.1) {
  Dico.norm <- CreationBases(data$position)
  
  t.cont <- seq(min(data$position), max(data$position), by = byseq)
  t.obs <- sort(unique(data$position))
  # browser()
  df.F <- data.frame(
    c(t.cont, t.cont),
    c(f.hat.old(
      t = t.cont,
      coef = model$Coef.Val, group = model$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    ) - mean(f.hat.old(
      t = t.obs,
      coef = model$Coef.Val, group = model$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    )), f.hat.old(
      t = t.cont,
      coef = model$Coef.Val, group = 1 - model$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    ) - mean(f.hat.old(
      t = t.obs,
      coef = model$Coef.Val, group = 1 - model$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    ))),
    c(rep(0, length(t.cont)), rep(1, length(t.cont)))
  )
  
  colnames(df.F) <- c("position", "F.fit", "Group")
  return(df.F)
}

create.CI <- function(diff.CI, data, sig.init) {
  dbar <- colMeans(do.call("rbind", lapply(diff.CI, function(x) {
    x$diff
  })))
  
  sbar <- sqrt(colMeans(do.call("rbind", lapply(diff.CI, function(x) {
    (x$diff - dbar)^2
  }))))
  
  Mb <- unlist(lapply(diff.CI, function(x) {
    max(abs(x$diff - dbar) / sbar)
  }))
  
  qb <- quantile(Mb, probs = 0.975)
  
  CIlow <- data.frame(position = diff.CI[[1]]$position, low = dbar - qb * sbar)
  CIup <- data.frame(position = diff.CI[[1]]$position, up = dbar + qb * sbar)
  
  set.seed(123)
  obs = diff.f(pred.f(fit.boot(data, sig.init), data, byseq = 0.1))
  
  df.f <- data.frame(obs, CIlow[, 2], CIup[, 2])
  colnames(df.f) <- c("Month", "Group diff.", "CI lower", "CI upper")
  # df.f = apply(df.f, 2, function(x) round(x, 2))
  rownames(df.f) <- NULL
  
  return(as.data.frame(df.f))
}

f.test <- function(x, y, series, t,  name_group_var = "group", plmm_output, n_boot = 1000) {
  # data.f <- data
  y <- y - plmm_output$lasso_output$x_fit - rep(plmm_output$out_phi$phi, plmm_output$ni)
  # data.f <- data.f[, c("series", "position", "Y", "Group")]
  
  t.obs <- sort(unique(t))
  
  data <- data.frame(y, series, t, x[,name_group_var])
  colnames(data)[4] = name_group_var
  
  samples <- boot_samples(data = data, n_boot = n_boot)
  
  pb <- utils::txtProgressBar(min = 0, max = length(samples), style = 3)
  
  pre.boot <- boot_pre(data)
  sig.init <-scalreg::scalreg(scale(pre.boot), scale(data.f$Y))$hsigma
  
  res.boot <- list()
  for (k in 1:length(samples)) {
    res.boot[[k]] <- fit.boot(Data = samples[[k]], sig.init)
    
    utils::setTxtProgressBar(pb, k)
  }
  
  cat("\nCompleted fitting Bootstrap samples. Now formatting results, and generating figure.\n")
  
  L2norm.df <- overall.test(res.boot, plmm_output = plmm_output)
  
  pred.CI <- lapply(1:length(res.boot), function(i) {
    pred.f(res.boot[[i]], samples[[i]])
  })
  
  
  diff.CI <- lapply(pred.CI, diff.f)
  df.CI <- create.CI(diff.CI, data.f, sig.init)
  
  plot.CI <- ggplot(df.CI, aes(x = Month, y = `Group diff.`)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(aes(Month, `CI lower`), data = df.CI, col = "blue") +
    geom_line(aes(Month, `CI upper`), data = df.CI, col = "blue") +
    labs(x = "Time", y = "Difference") +
    scale_x_continuous(breaks = t.obs)
  
  return(list(L2norm.df = L2norm.df, df.CI = df.CI, plot.CI = plot.CI))
}
