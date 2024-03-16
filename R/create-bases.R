filter_nonzero_bases <- function(bases) {
  # Filter bases functions that are essentially 0 if any
  bases.Present <- c()
  M <- ncol(bases)

  PresentBases <- c()
  for (j in 1:M) {
    fj <- bases[, j]
    if (sum(abs(fj)) > 10^-10) {
      PresentBases <- c(PresentBases, j)
    }
  }
  bases.Present <- bases[, PresentBases]
  return(list(Fh.P = bases.Present, PresentBases = PresentBases))
}

create_bases <- function(position, keep = NULL) {
  # Inputs:
  # position: time for all individuals
  #
  # Output: a list with:
  # - F.Bases: Functions of the dictionary
  # - Num.Bases.Pres: ID of the functions

  # Choosing the dictionary for Lasso estimation

  # Fourier (sin,cos)
  Ff <- c()

  n <- max(position)
  Ff <- c()
  for (i in seq(1, n, 1)) {
    Ff <- cbind(Ff, sin(2 * pi * i * position / n), cos(2 * pi * i * position / n))
  }

  # Power functions
  x <- position / n
  Fx <- c()
  rg <- seq(0.1, 2, 0.02)
  for (i in 1:length(rg)) {
    Fx <- cbind(Fx, x^rg[i])
  }

  F.tot <- cbind(Ff, Fx)

  if (!is.null(keep)) {
    F.tot <- F.tot[, keep]
  }

  F.Bases <- c()
  F.Bases.tot <- c()
  F.Bases.tot <- filter_nonzero_bases(F.tot)
  Num.Bases.Pres <- F.Bases.tot$PresentBases
  F.Bases <- F.Bases.tot$Fh.P

  return(list(F.Bases = F.Bases, Num.Bases.Pres = Num.Bases.Pres))
}
