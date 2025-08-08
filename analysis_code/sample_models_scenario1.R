library(dplyr)
library(purrr)


devtools::load_all()

models <- load_FvCBmodels()
# Grids
Ci     <- c(20,50,100,150,200,300,400,500,600,700,800,1000,1200)
PPFDs  <- c(20,100,300,500,700,1000,1200,1500)

envs_template <- list(
  C_i = Ci,
  O   = rep(210, length(Ci))
)

# Mean parameters # please modify it accordingly
param_mean <- c(
  Vcmax = 70,
  Jmax  = 130,
  J     = 130,
  Rd    = 2,
  KC    = 270,
  KO    = 165,
  Gstar = 38,
  alphaJ = 0.3,
  thetaJ = 0.8,
  Vtpu   = 9,
  gm     = 0.3,
  atpu   = 0.3,
  phi2m   = 0.85,  # will map to phi2m where needed
  fQ     = 1,
  fpseudo = 0,
  fcyc    = NA_real_, # if NA we set to 1 - fQ - fpseudo
  h       = 14
)

sample_params <- function(param_mean, pct_range = 0.1) {
  p <- sapply(param_mean, function(mu) {
    lo <- mu * (1 - pct_range)
    hi <- mu * (1 + pct_range)
    runif(1, min = lo, max = hi)
  })

  p["fcyc"] <- NA_real_

  return(p)
}


# Run N samples for ONE model -> 3D array [Ci, PPFD, samples]
simulate_model_samples <- function(model_def, N, envs_template, Ci, PPFDs,param_mean){
  n_ci   <- length(Ci)
  n_ppfd <- length(PPFDs)
  out <- array(NA_real_, dim = c(n_ci, n_ppfd, N),
               dimnames = list(Ci = as.character(Ci),
                               PPFD = as.character(PPFDs),
                               sample = paste0("s", seq_len(N))))
  for (s in seq_len(N)) {
    p <- sample_params(param_mean, 0.1)
    pars <- model_def$map(p)

    for (j in seq_along(PPFDs)) {
      envs <- envs_template
      envs$PPFD <- rep(PPFDs[j], length(Ci))

      res <- model_def$fun(envs = envs, pars = pars)
      A   <- res$An
      if (length(A) != n_ci) {
        stop("Model returned A of length ", length(A), " but expected ", n_ci)
      }
      out[, j, s] <- A
    }
  }
  return(out)
}

# --------- RUN IT ----------
N <- 10000  # number of parameter samples

result_arrays <- imap(models, ~ simulate_model_samples(.x, N, envs_template, Ci, PPFDs,param_mean))
# result_arrays is a named list; each element is a 3D array [Ci, PPFD, N].
# Example access: result_arrays[["FvCB (1980)"]][ , "500", 1]
