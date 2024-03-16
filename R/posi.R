library(dplyr)

debias.plmm <- function(simu, model, a = 1, level = 0.95, Z = NULL) {
  data.lmm <- simu
  data.lmm$Y <- data.lmm$Y - model$Res.F$out.F$F.fit
  
  X <- as.matrix(subset(data.lmm, select = -c(series, position, Y)))
  y <- data.lmm$Y
  
  grp <- factor(data.lmm$series)
  
  z <- model.matrix(~ as.factor(series) - 1, data.lmm[, "series", drop = F])
  
  Sigma <- a * model$su * z %*% t(z) + model$se * diag(rep(1, length(y)))
  
  Sigma.svd <- svd(Sigma)
  Sig.a.inv.half <- Sigma.svd$u %*% diag(1 / sqrt(Sigma.svd$d)) %*% t(Sigma.svd$u)
  
  X.a <- Sig.a.inv.half %*% X
  y.a <- Sig.a.inv.half %*% y
  
  if(is.null(Z)) {
    X.a <- X.a[, colnames(X.a) %in% c(
      names(model$Res.F$theta[-1][model$Res.F$theta[-1] != 0]),
      sample(names(model$Res.F$theta[-1][model$Res.F$theta[-1] == 0]), 5)
    )]
    
    de.sparsified.sl <- hdi::lasso.proj(X.a, y.a,
                                        suppress.grouptesting = TRUE, return.Z = T,
                                        do.ZnZ = T, betainit = "scaled lasso"
    )
    
    lp.Z <- de.sparsified.sl$Z
  } else {
    lp.Z = Z
  }
  
  
  N <- length(y)
  n <- length(unique(grp))
  q <- ncol(z)
  p <- ncol(X)
  
  beta.db.mlm <- NULL
  beta.db.sd.mlm <- NULL
  
  beta.hat <- model$Res.F$theta[names(model$Res.F$theta) %in% colnames(X.a)]
  
  res <- y.a - X.a %*% beta.hat
  
  for (j in 1:length(beta.hat)) {
    col.j <- j
    
    wj.mlm <- lp.Z[, j]
    
    beta.db.mlm[j] <- beta.hat[col.j] + sum(wj.mlm * res) / sum(wj.mlm * X.a[, col.j])
    
    ok <- data.frame(wj.mlm, res, grp) %>%
      group_by(grp) %>%
      summarize(num = (sum(wj.mlm * res))^2)
    num <- sum(ok$num)
    
    beta.db.sd.mlm[j] <- sqrt(num) / sqrt((sum(wj.mlm * X.a[, col.j]))^2)
  }
  
  if (level == 0.95) {
    ci <- cbind(
      beta.db.mlm - 1.96 * beta.db.sd.mlm,
      beta.db.mlm + 1.96 * beta.db.sd.mlm
    )
  }
  
  if (level == 0.90) {
    ci <- cbind(
      beta.db.mlm - 1.64 * beta.db.sd.mlm,
      beta.db.mlm + 1.64 * beta.db.sd.mlm
    )
  }
  
  if (level == 0.99) {
    ci <- cbind(
      beta.db.mlm - 2.58 * beta.db.sd.mlm,
      beta.db.mlm + 2.58 * beta.db.sd.mlm
    )
  }
  
  # Calculate p-value
  bprojrescaled <- beta.db.mlm * (1 / beta.db.sd.mlm)
  pv <- 2 * pnorm(abs(bprojrescaled), lower.tail = FALSE)
  
  names(beta.db.mlm) <- colnames(X.a)
  
  debias.PLMM <- list(
    bhat = beta.db.mlm, bhat.sd = beta.db.sd.mlm, pv = pv,
    ci = ci, beta.hat = beta.hat
  )
  
  df.posi = cbind(debias.PLMM$beta.hat, debias.PLMM$bhat, debias.PLMM$ci, debias.PLMM$pv)
  colnames(df.posi) = c("Estimate", "Debiased", "CI lower", "CI upper", "p-value")
  df.posi[, 1:4] <- round(df.posi[, 1:4], 2)
  df.posi = as.data.frame(df.posi)
  df.posi[,5] = roundpval(df.posi[,5])
  
  return(df.posi)
}
