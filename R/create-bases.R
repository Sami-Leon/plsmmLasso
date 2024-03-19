filter_nonzero_bases <- function(bases) {
  # Filter bases functions that are essentially 0 if any
  filtered_bases <- c()
  M <- ncol(bases)
  bases_indices <- c()
  for (j in 1:M) {
    fj <- bases[, j]
    if (sum(abs(fj)) > 10^-10) {
      bases_indices <- c(bases_indices, j)
    }
  }
  filtered_bases <- bases[, bases_indices]
  return(list(filtered_bases = filtered_bases, selected_bases = bases_indices))
}

create_bases <- function(t, keep = NULL) {
  max_t <- max(t)
  n_timepoints <- length(t)
  
  Fourier_basis <- matrix(NA, n_timepoints, 2 * max_t)
  for (i in 1:max_t) {
    angle <- 2 * pi * i * t / max_t
    Fourier_basis[, (2 * i - 1):(2 * i)] <- cbind(sin(angle), cos(angle))
  }

  normalized_t <- t / max_t
  poly_degrees <- seq(0.1, 2, 0.02)
  
  poly_basis <- matrix(NA, n_timepoints, length(poly_degrees))
  for (i in 1:length(poly_degrees)) {
    poly_basis[,i] <- normalized_t^poly_degrees[i]
  }

  bases_functions <- cbind(Fourier_basis, poly_basis)

  if (!is.null(keep)) {
    bases_functions <- bases_functions[, keep]
  }

  filtered_bases_functions <- filter_nonzero_bases(bases = bases_functions)
  
  selected_bases <- filtered_bases_functions$selected_bases
  filtered_bases_bases_functions <- filtered_bases_functions$filtered_bases

  return(list(bases = filtered_bases_bases_functions, selected_bases = selected_bases))
}
