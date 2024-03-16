filter_nonzero_bases <- function(bases) {
  # Filter bases functions that are essentially 0 if any
  filted_bases <- c()
  M <- ncol(bases)

  bases_indices <- c()
  for (j in 1:M) {
    fj <- bases[, j]
    if (sum(abs(fj)) > 10^-10) {
      bases_indices <- c(bases_indices, j)
    }
  }
  filted_bases <- bases[, bases_indices]
  return(list(Fh.P = filted_bases, PresentBases = bases_indices))
}

create_bases <- function(t, keep = NULL) {
  max_t <- max(t)
  n_timepoints <- length(t)
  
  Fourier_basis <- matrix(NA, n_timepoints, 2 * max_t)
  for (i in 1:max_t) {
    angle <- 2 * pi * i * t / max_t
    Fourier_basis <- cbind(Fourier_basis, sin(angle), cos(angle))
  }

  normalized_t <- t / max_t
  poly_degrees <- seq(0.1, 2, 0.02)
  
  poly_basis <- matrix(NA, n_timepoints, length(poly_degrees))
  for (i in 1:length(poly_degrees)) {
    poly_basis <- cbind(poly_basis, normalized_t^poly_degrees[i])
  }

  bases_functions <- cbind(Fourier_basis, poly_basis)

  if (!is.null(keep)) {
    bases_functions <- bases_functions[, keep]
  }

  filtered_bases_functions <- filter_nonzero_bases(bases_functions)
  
  bases_indices <- filtered_bases_functions$PresentBases
  bases_functions_out <- filtered_bases_functions$Fh.P

  return(list(F.Bases = bases_functions_out, Num.Bases.Pres = bases_indices))
}
