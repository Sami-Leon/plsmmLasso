sample_boot <- function(data, n_boot) {
  unique_series <- unique(data$series)
  n_series <- length(unique_series)

  bootstrap_samples <- lapply(1:n_boot, function(x) {
    sample_id <- sample(unique_series, n_series, replace = TRUE)
    mat_boot <- do.call(rbind, lapply(sample_id, function(j) data[data$series == j, ]))
    return(mat_boot)
  })

  return(bootstrap_samples)
}

create_bases_boot <- function(data, timexgroup = TRUE) {
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
    bases_timexgroup[!index_ref_group, (M + 1):(2 * M)] <- bases[!index_ref_group, ]

    bases <- bases_timexgroup
  }


  return(bases)
}

fit_boot <- function(data, min_lambda) {
  boot_bases <- create_bases_boot(data)

  lambda_grid <- seq(1.2, min_lambda, -0.1) * sqrt(2 * log(ncol(boot_bases)) / length(data$y))

  cv_fit <- glmnet::cv.glmnet(boot_bases, data$y,
    alpha = 1,
    lambda = lambda_grid
  )

  final_fit <- glmnet::glmnet(boot_bases, data$y, alpha = 1, lambda = cv_fit$lambda.min)
  fitted_values <- glmnet::predict.glmnet(final_fit, newx = boot_bases, s = cv_fit$lambda.min)

  data$f_fit <- fitted_values

  out_f <- data[, c("t", "f_fit", "group")]

  out_f[out_f$Group == 0, ]$f_fit <- out_f[out_f$Group == 0, ]$f_fit - attr(scale(unique(out_f[out_f$Group == 0, ]$f_fit),
    scale = FALSE
  ), "scaled:center")
  out_f[out_f$Group == 1, ]$f_fit <- out_f[out_f$Group == 1, ]$f_fit - attr(scale(unique(out_f[out_f$Group == 1, ]$f_fit),
    scale = FALSE
  ), "scaled:center")
  return(list(out_f = out_f, alpha = as.vector(final_fit$beta[, 1])))
}

calc_f_diff <- function(out_f) {
  fit0 <- out_f[out_f$group == 0, ]
  fit0 <- fit0[!duplicated(fit0$t), ]
  fit0 <- fit0[order(fit0$t), ]

  fit1 <- out_f[out_f$group == 1, ]
  fit1 <- fit1[!duplicated(fit1$t), ]
  fit1 <- fit1[order(fit1$t), ]

  out <- data.frame(t = fit1$t, diff = fit0$f_fit - fit1$f_fit)
  return(out)
}

test_overall_f <- function(list_fitted_boot, plmm_output) {
  df_list_fit <- lapply(seq_along(list_fitted_boot), function(i) {
    df <- list_fitted_boot[[i]]$out_f
    df$boot <- as.character(i)
    return(df)
  })

  diff_list <- lapply(df_list_fit, calc_f_diff)

  T_obs <- sum((calc_f_diff(plmm_output$lasso_output$out_f)$diff)^2)
  T_boot <- sapply(diff_list, function(x) {
    sum(x$diff^2)
  })

  pvalue <- pnorm(abs((2 * T_obs - mean(T_boot)) / sqrt(var(T_boot))),
    lower.tail = FALSE
  ) * 2

  overall_f <- data.frame(T_obs, pvalue)

  colnames(overall_f) <- c("T", "p-value")

  return(overall_f)
}

pred_f <- function(model, data, byseq = 0.1) {
  selected_bases <- create_bases(data$t)$selected_bases

  t_cont <- seq(min(data$t), max(data$t), by = byseq)
  t_obs <- sort(unique(data$t))

  df.F <- data.frame(
    c(t_cont, t_cont),
    c(f_predict(
      t = t_cont,
      coef = model$alpha, group = model$out_f$group[1],
      keep = selected_bases
    ) - mean(f_predict(
      t = t_obs,
      coef = model$alpha, group = model$out_f$group[1],
      keep = selected_bases
    )), f_predict(
      t = t_cont,
      coef = model$alpha, group = 1 - model$out_f$group[1],
      keep = selected_bases
    ) - mean(f_predict(
      t = t_obs,
      coef = model$alpha, group = 1 - model$out_f$group[1],
      keep = selected_bases
    ))),
    c(rep(0, length(t_cont)), rep(1, length(t_cont)))
  )

  colnames(df.F) <- c("t", "f_fit", "group")
  return(df.F)
}

create_CI <- function(list_diff_CI, data, min_lambda) {
  d_bar <- colMeans(do.call("rbind", lapply(list_diff_CI, function(x) {
    x$diff
  })))

  s_bar <- sqrt(colMeans(do.call("rbind", lapply(list_diff_CI, function(x) {
    (x$diff - d_bar)^2
  }))))

  M_b <- unlist(lapply(list_diff_CI, function(x) {
    max(abs(x$diff - d_bar) / s_bar)
  }))

  q_b <- quantile(M_b, probs = 0.975)

  CI_low <- data.frame(t = list_diff_CI[[1]]$t, low = d_bar - q_b * s_bar)
  CI_up <- data.frame(t = list_diff_CI[[1]]$t, up = d_bar + q_b * s_bar)

  obs <- calc_f_diff(pred_f(fit_boot(data, min_lambda), data, byseq = 0.1))

  CI_diff_f <- data.frame(obs, CI_low[, 2], CI_up[, 2])
  colnames(CI_diff_f) <- c("t", "f diff.", "Lower 95%", "Upper 95%")
  rownames(CI_diff_f) <- NULL

  return(CI_diff_f)
}

test_f <- function(x, y, series, t, name_group_var, plmm_output, n_boot = 1000,
                   predicted = FALSE) {
  y <- y - plmm_output$lasso_output$x_fit - rep(plmm_output$out_phi$phi, plmm_output$ni)

  t_obs <- sort(unique(t))

  data <- data.frame(y, series, t, x[, name_group_var])
  colnames(data)[4] <- "group"

  samples <- sample_boot(data = data, n_boot = n_boot)

  pb <- utils::txtProgressBar(min = 0, max = length(samples), style = 3)

  min_lambda <- scalreg::scalreg(scale(create_bases_boot(data)), scale(data$y))$hsigma

  fitted_boot <- vector("list", n_boot)

  for (k in 1:n_boot) {
    fitted_boot[[k]] <- fit_boot(data = samples[[k]], min_lambda = min_lambda)

    utils::setTxtProgressBar(pb, k)
  }

  cat("\nCompleted fitting Bootstrap samples. Now formatting results, and generating figure.\n")

  overall_test_results <- test_overall_f(
    list_fitted_boot = fitted_boot,
    plmm_output = plmm_output
  )

  predicted_f <- lapply(1:length(fitted_boot), function(i) {
    pred_f(model = fitted_boot[[i]], data = samples[[i]])
  })


  diff_predicted_f <- lapply(predicted_f, calc_f_diff)
  CI_f <- create_CI(diff_predicted_f, data, min_lambda)

  if (predicted) {
    plot_CI <- ggplot2::ggplot(CI_f, aes(x = t, y = `f diff.`)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_ribbon(aes(x = t, ymin = `Lower 95%`, ymax = `Upper 95%`),
        data = CI_f,
        fill = "gray", alpha = 0.6
      ) +
      geom_line(linewidth = 0.7) +
      geom_point(
        aes(x = t, y = `f diff.`),
        CI_f[CI_f$t %in% t_obs, ]
      ) +
      scale_x_continuous(breaks = t_obs)
  } else {
    CI_obs_f <- CI_f[CI_f$t %in% t_obs, ]

    plot_CI <- ggplot2::ggplot(CI_obs_f, aes(x = t, y = `f diff.`)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_ribbon(aes(x = t, ymin = `Lower 95%`, ymax = `Upper 95%`),
        data = CI_obs_f,
        fill = "gray", alpha = 0.6
      ) +
      geom_line(linewidth = 0.7) +
      geom_point() +
      scale_x_continuous(breaks = t_obs)
  }

  return(list(
    overall_test_results = overall_test_results, CI_f = CI_f,
    plot_CI = plot_CI
  ))
}