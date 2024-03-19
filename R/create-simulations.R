f <- function(t, A, omega, cos = FALSE) {
  if (cos) {
    return(as.vector(-A * sin(2 * pi * t / omega)))
  }
  as.vector(A * sin(2 * pi * t / omega))
}

cst_cor <- function(n, rho) {
  m <- matrix(rho, n, n)
  diag(m) <- 1
  return(m)
}

simulate_group_inter <- function(N = 50, n_mvnorm = 100, grouped = TRUE,
                                 timepoints = 3:5, nonpara_inter = TRUE,
                                 sample_from, cst_ni, cos = FALSE, A_vec = c(1, 1.5)) {
  if (nonpara_inter) {
    A <- c(A_vec[1], A_vec[2])
    omega <- c(60, 110)
    
    f0_mean <- mean(f(sample_from, A[1], omega[1], cos = FALSE))
    f1_mean <- mean(f(sample_from, A[2], omega[2], cos = cos))
  } else {
    A <- c(A_vec[1], A_vec[1])
    omega <- c(60, 60)
    f0_mean <- mean(f(sample_from, A[1], omega[1], cos = FALSE))
    f1_mean <- f0_mean
  }
  
  phi <- rnorm(N, 0, sqrt(0.5))
  
  y <- NULL
  sim <- NULL
  f_val <- NULL
  
  for (i in 1:N) {
    if (cst_ni) {
      ni <- timepoints
    } else {
      ni <- sample(timepoints, 1)
    }
    
    if (grouped) {
      theta <- c(3, 2, 1)
    } else {
      theta <- c(0, 2, 0)
    }
    
    group <- rep(sample(c(0, 1), 1), ni)
    
    x1 <- rep(rnorm(1, 1, sqrt(0.5)), ni)
    
    eps <- rnorm(ni, 0, sqrt(0.2))
    
    t <- sort(sample(sample_from, ni, replace = F))
    
    if (group[1] == 0) {
      sim <- rbind(sim, cbind(
        rep(i, ni), t, phi[i] + f(t, A[1], omega[1], cos = FALSE) + eps - f0_mean,
        group, x1
      ))
      
      f_val <- c(f_val, f(t, A[1], omega[1], cos = FALSE) - f0_mean)
    } else {
      sim <- rbind(sim, cbind(
        rep(i, ni), t, phi[i] + f(t, A[2], omega[2], cos = cos) + eps - f1_mean,
        group, x1
      ))
      
      f_val <- c(f_val, f(t, A[2], omega[2], cos = cos) - f1_mean)
    }
  }
  
  x <- MASS::mvrnorm(nrow(sim), rep(0, n_mvnorm - 1), cst_cor(n_mvnorm - 1, 0))
  
  sim <- cbind(sim, x)
  
  colnames(sim) <- c("series", "t", "y", "group", paste0("x", 1:(ncol(x) + 1)))
  
  sim[, "x2"] <- sim[, "group"] * sim[, "x1"]
  sim[, "y"] <- sim[, "y"] + sim[, c("group", "x1", "x2"), drop = F] %*% theta
  
  sim <- as.data.frame(sim)
  
  phi <- rep(phi, table(sim$series))
  
  sim <- sim[order(sim$series, sim$t), ]
  
  f_val <- f_val[order(sim$series, sim$t)]
  
  phi <- phi[order(sim$series, sim$t)]
  
  return(list(sim = sim, phi = phi, f_val = f_val))
}

f_predict = function(t, coef, group, keep = NULL) {
  
  bases = create_bases(t, keep = keep)$bases
  
  if(group == 0) {
    coef = coef[1:ncol(bases)]
  } else {
    coef = coef[(ncol(bases)+1):length(coef)]
  }
  
  return(bases %*% coef)
  
}

plot_fit <- function(x, y, series, t,  name_group_var = "group", 
                     plmm_output, predicted = FALSE) {
  data <- data.frame(y, series, t, x)
  
  data$f_fit <- plmm_output$lasso_output$out_f$f_fit
  data$x_fit <- plmm_output$lasso_output$x_fit
  data$phi <- rep(plmm_output$out_phi$phi, table(data$series))

  bases_functions <- create_bases(t)

  t_obs = sort(unique(t))
  
  t_cont <- seq(min(t_obs), max(t_obs), by = 0.1)
  
  predicted_f <- data.frame(c(t_cont, t_cont),
    c(f_predict(
      t = t_cont,
      coef = plmm_output$lasso_output$alpha, group = plmm_output$lasso_output$out_f$group[1],
      keep = bases_functions$selected_bases
    ) - mean(f_predict(
      t = t_obs,
      coef = plmm_output$lasso_output$alpha, group = plmm_output$lasso_output$out_f$group[1],
      keep = bases_functions$selected_bases
    )), f_predict(
      t = t_cont,
      coef = plmm_output$lasso_output$alpha, group = 1 - plmm_output$lasso_output$out_f$group[1],
      keep = bases_functions$selected_bases
    ) - mean(f_predict(
      t = t_obs,
      coef = plmm_output$lasso_output$alpha, group = 1 - plmm_output$lasso_output$out_f$group[1],
      keep = bases_functions$selected_bases
    ))),
    c(rep(0, length(t_cont)), rep(1, length(t_cont)))
  )
  
  colnames(predicted_f) <- c("t", "f_cont", "group")
  
  means <- aggregate(cbind(phi, x_fit) ~ group, data = data, FUN = mean)
  names(means) <- c("group", "phi", "x_fit")
  
  predicted_f <- merge(predicted_f, means, by = "group")
  
  predicted_f$mean_trajectories <- predicted_f$f_cont + predicted_f$x_fit + predicted_f$phi
  
  obs_f = predicted_f[predicted_f$t %in% t_obs, ]
  
  p <- ggplot2::ggplot(data = data, aes(x = t, y = y))
  
  if(predicted) {
    p.F.overall <- p + geom_line(aes(x = t, y = y, group = series)) +
      geom_line(aes(x = t, y = mean_trajectories), data = predicted_f, size = 1,
                col = "red") +
      facet_grid(. ~ group) + geom_point(aes(x = t, y = mean_trajectories), 
                                         data = obs_f, size = 2,
                                         col = "red") +
      scale_x_continuous(breaks = t_obs)
    
    p.F <- ggplot(aes(x = t, y = f_cont), data = predicted_f) +
      geom_line(size = 1, col = "red") +
      facet_grid(. ~ group) +
      geom_point(aes(x = t, y = f_cont),
                 data = obs_f, size = 2, col = "red") +
      scale_x_continuous(breaks = t_obs)
  } else {
    p.F.overall <- p + geom_line(aes(x = t, y = y, group = series)) +
      geom_line(aes(x = t, y = mean_trajectories), data = obs_f, size = 1,
                col = "red") +
      geom_point(aes(x = t, y = mean_trajectories), 
                 data = obs_f, size = 2, col = "red") +
      facet_grid(. ~ group) + 
      scale_x_continuous(breaks = t_obs)
    
    p.F <- ggplot(aes(x = t, y = f_cont), data = obs_f) +
      geom_line(size = 1, col = "red") +
      facet_grid(. ~ group) +
      geom_point(size = 2, col = "red") +
      scale_x_continuous(breaks = t_obs)
  }
  
  ggarrange(p.F.overall, p.F, ncol = 1, nrow = 2, common.legend = TRUE, 
            legend = "bottom")
}
