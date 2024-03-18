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
  
  colnames(sim) <- c("series", "position", "y", "group", paste0("x", 1:(ncol(x) + 1)))
  
  sim[, "x2"] <- sim[, "group"] * sim[, "x1"]
  sim[, "y"] <- sim[, "y"] + sim[, c("group", "x1", "x2"), drop = F] %*% theta
  
  sim <- as.data.frame(sim)
  
  phi <- rep(phi, table(sim$series))
  
  sim <- sim[order(sim$series, sim$position), ]
  
  f_val <- f_val[order(sim$series, sim$position)]
  
  phi <- phi[order(sim$series, sim$position)]
  
  return(list(sim = sim, phi = phi, f_val = f_val))
}

f.hat.old = function(t, coef, group, keep = NULL) {
  
  F.Bases = CreationBases(t, keep = keep)$F.Bases
  
  if(group == 0) {
    coef = coef[1:ncol(F.Bases)]
  } else {
    coef = coef[(ncol(F.Bases)+1):length(coef)]
  }
  
  return(F.Bases %*% coef)
  
}

plot.fit <- function(list.EM.out, data, same = FALSE) {
  EM.out <- list.EM.out
  
  data1 <- data[[1]]
  data1$phi <- data[[2]]
  data1$f <- data[[3]]
  
  t.obs <- sort(unique(data1$position))
  
  data1$X <- as.matrix(data1[, 4:6]) %*% cbind(c(3, 2, 1))
  
  bluemean_X <- data1 %>%
    group_by(Group) %>%
    mutate(mean = mean(X)) %>%
    ungroup()
  
  bluemean_U2 <- data1 %>%
    group_by(Group) %>%
    mutate(mean = mean(phi)) %>%
    ungroup()
  
  data1$blueline <- data1$f + bluemean_U2$mean + bluemean_X$mean
  
  mean_blueline <- data1 %>%
    group_by(Group, position) %>%
    mutate(mean = mean(blueline)) %>%
    ungroup()
  
  data1$blueline <- mean_blueline$mean
  
  data1$F.fit <- EM.out$Res.F$out.F$F.fit
  data1$X.fit <- EM.out$Res.F$X.fit
  data1$U2 <- rep(EM.out$U2$U2, table(data1$series))
  
  mean_X <- data1 %>%
    group_by(Group) %>%
    mutate(mean = mean(X.fit)) %>%
    ungroup()
  
  mean_U2 <- data1 %>%
    group_by(Group) %>%
    mutate(mean = mean(U2)) %>%
    ungroup()
  
  redline <- EM.out$Res.F$out.F$F.fit + mean_X$mean + mean_U2$mean
  data1$redline <- redline
  
  mean_redline <- data1 %>%
    group_by(Group, position) %>%
    mutate(mean = mean(redline)) %>%
    ungroup()
  
  data1$redline <- mean_redline$mean
  
  p <- ggplot(data = data1, aes(x = position, y = Y))
  
  group_label <- c("0" = "Group1", "1" = "Group2")
  
  F.plot <- p + geom_line(aes(x = position, y = Y, group = series)) +
    geom_line(aes(x = position, y = blueline, color = "truth"), data = data1, size = 1.3) +
    geom_line(aes(x = position, y = redline, color = "estimate"),
              data = data1, size = 1.3
    ) +
    scale_color_manual(name = "Legend", values = c("truth" = "blue", "estimate" = "red")) +
    facet_grid(. ~ Group, labeller = as_labeller(group_label))
  
  p <- ggplot(data = data1, aes(x = position, y = Y))
  
  group_label <- c("0" = "Group1", "1" = "Group2")
  
  data1.dedup <- data1[!duplicated(data1$series), ]
  data1.dedup <- data1.dedup[, c("Group", "phi", "U2")]
  colnames(data1.dedup)[2:3] <- c("truth", "estimate")
  
  Dico.norm <- CreationBases(data1$position)
  
  
  t.cont <- seq(min(t.obs), max(t.obs), by = 0.1)
  values <- f(t.cont, A = 1, omega = 60) - mean(f(t.obs, A = 1, omega = 60))
  
  if (same) {
    values.1 <- values
  } else {
    values.1 <- f(t.cont, A = 1.5, omega = 110) - mean(f(t.obs, A = 1.5, omega = 110))
  }
  
  
  
  df.F <- data.frame(
    c(t.cont, t.cont), c(values, values.1),
    c(f.hat.old(
      t = t.cont,
      coef = EM.out$Res.F$Coef.Val, group = EM.out$Res.F$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    ) - mean(f.hat.old(
      t = t.obs,
      coef = EM.out$Res.F$Coef.Val, group = EM.out$Res.F$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    )), f.hat.old(
      t = t.cont,
      coef = EM.out$Res.F$Coef.Val, group = 1 - EM.out$Res.F$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    ) - mean(f.hat.old(
      t = t.obs,
      coef = EM.out$Res.F$Coef.Val, group = 1 - EM.out$Res.F$out.F$Group[1],
      keep = Dico.norm$Num.Bases.Pres
    ))),
    c(rep(0, length(t.cont)), rep(1, length(t.cont)))
  )
  
  colnames(df.F) <- c("t", "f", "F.fit", "Group")
  
  means <- aggregate(cbind(phi, X, U2, X.fit) ~ Group, data = data1, FUN = mean)
  names(means) <- c("Group", "phi", "X", "U2", "X.fit")
  
  df.F <- merge(df.F, means, by = "Group")
  
  df.F$f.overall <- df.F$f + df.F$X + df.F$phi
  df.F$F.fit.overall <- df.F$F.fit + df.F$X.fit + df.F$U2
  
  
  p.F.overall <- p + geom_line(aes(x = position, y = Y, group = series)) +
    geom_line(aes(x = t, y = f.overall, color = "truth"), size = 1.3, data = df.F) +
    geom_line(aes(x = t, y = F.fit.overall, color = "estimate"), data = df.F, size = 1.3) +
    scale_color_manual(name = "Legend", values = c("truth" = "blue", "estimate" = "red")) +
    facet_grid(. ~ Group, labeller = as_labeller(group_label)) +
    ylab("Y") + xlab("Time") + geom_point(aes(x = t, y = f.overall), col = "blue", data = df.F[df.F$t %in% t.obs, ], size = 2.5) +
    geom_point(aes(x = t, y = F.fit.overall), col = "red", data = df.F[df.F$t %in% t.obs, ], size = 2.5) +
    scale_x_continuous(breaks = t.obs)
  
  p.F <- ggplot(aes(x = t, y = f), data = df.F) +
    geom_line(aes(x = t, y = f, color = "truth"), size = 1.3, data = df.F) +
    geom_line(aes(x = t, y = F.fit, color = "estimate"), data = df.F, size = 1.3) +
    scale_color_manual(name = "Legend", values = c("truth" = "blue", "estimate" = "red")) +
    facet_grid(. ~ Group, labeller = as_labeller(group_label)) +
    ylab("Y") +
    xlab("Time") +
    geom_point(aes(x = t, y = f), col = "blue", data = df.F[df.F$t %in% t.obs, ], size = 2.5) +
    geom_point(aes(x = t, y = F.fit), col = "red", data = df.F[df.F$t %in% t.obs, ], size = 2.5) +
    ylim(c(-1.5, 10)) +
    scale_x_continuous(breaks = t.obs)
  
  ggarrange(p.F.overall, p.F, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
}
