debias_plmm <- function(x, y, series, plmm_output, a = 1, Z = NULL) {
  y_offset <- y - plmm_output$lasso_output$out_f$f_fit

  series <- as.factor(series)

  z <- model.matrix(~ series - 1)

  Sigma_a <- a * plmm_output$su * z %*% t(z) + plmm_output$se * diag(rep(1, length(y_offset)))

  Sigma_a_svd <- svd(Sigma_a)
  Sigma_a_sqrt_inv <- Sigma_a_svd$u %*% diag(1 / sqrt(Sigma_a_svd$d)) %*% t(Sigma_a_svd$u)

  x_a <- Sigma_a_sqrt_inv %*% x
  y_a <- Sigma_a_sqrt_inv %*% y_offset

  if (is.null(Z)) {
    de_sparsified <- hdi::lasso.proj(x_a, y_a,
      suppress.grouptesting = TRUE, return.Z = TRUE,
      do.ZnZ = TRUE, betainit = "scaled lasso"
    )

    debias_score_matrix <- de_sparsified$Z
  } else {
    debias_score_matrix <- Z
  }
  # removing intercept
  beta_original <- plmm_output$lasso_output$theta[-1] 

  res <- y_a - x_a %*% beta_original

  p <- length(beta_original)

  beta_debias <- rep(NA, p)
  beta_debias_sd <- rep(NA, p)

  for (j in 1:p) {
    score_j <- debias_score_matrix[, j]

    beta_debias[[j]] <- beta_original[j] + sum(score_j * res) / sum(score_j * x_a[, j])

    # Group by 'series' and summarize scaled_res
    df_res <- dplyr::summarize(dplyr::group_by(data.frame(score_j, res, series), series),
      scaled_res = (sum(score_j * res))^2
    )

    # Calculate scaled_rss
    scaled_rss <- sum(df_res$scaled_res)

    beta_debias_sd[[j]] <- sqrt(scaled_rss) / sqrt((sum(score_j * x_a[, j]))^2)
  }

  ci <- cbind(
    beta_debias - 1.96 * beta_debias_sd,
    beta_debias + 1.96 * beta_debias_sd
  )

  # Calculate p-value
  beta_debias_rescaled <- beta_debias * (1 / beta_debias_sd)
  pv <- 2 * pnorm(abs(beta_debias_rescaled), lower.tail = FALSE)

  df.posi <- data.frame(
    Estimate = beta_original,
    Debiased = beta_debias,
    `Std. Error` = beta_debias_sd,
    `Lower 95%` = ci[, 1],
    `Upper 95%` = ci[, 2],
    `p-value` = pv, check.names = FALSE
  )

  return(df.posi)
}
