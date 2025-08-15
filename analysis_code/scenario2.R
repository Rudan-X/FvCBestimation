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

# Define unique param_mean for each model 
# (Please change the values of the parameters accordingly to the values provided in the literature)
param_means <- list(
  "FvCB (1980)" = c(
    Vcmax = 70, Jmax = 130,  Rd = 2, KC = 270, KO = 165,
    Gstar = 38, alphaJ = 0.3, thetaJ = 0.8, Vtpu = 9
  ),
  "Harley (1992)" = c(
    Vcmax = 75, Jmax = 135,  Rd = 2.2, KC = 280, KO = 170,
    Gstar = 40,  gm = 0.32
  ),
  "Ethier & Long (2004)" = c(
    Vcmax = 68, Jmax = 128,  Rd = 1.8, KC = 265, KO = 160,
    Gstar = 37, alphaJ = 0.29, thetaJ = 0.79,  gm = 0.29
  ),
  "von Caemmerer (2000)" = c(
    Vcmax = 72, Jmax = 132,  Rd = 2.1, KC = 275, KO = 168,
    Gstar = 39, alphaJ = 0.31, thetaJ = 0.81, Vtpu = 9.2, atpu = 0.31
  ),
  "Yin (2004)" = c(
    Vcmax = 71, Jmax = 131,  Rd = 2.0, KC = 272, KO = 166,
    Gstar = 38.5,  thetaJ = 0.805,  phi2m = 0.855, 
    fQ = 1, fpseudo = 0, fcyc = NA_real_, h = 14.2
  ),
  "Dubois (2007)" = c(
    Vcmax = 73, J = 133, Rd = 2.3, KC = 278, KO = 172,
    Gstar = 41
  )
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

# Modified simulate_model_samples to accept param_mean as argument
simulate_model_samples <- function(model_def, N, envs_template, Ci, PPFDs, param_mean) {
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

# Modified to pass model-specific param_mean
result_arrays <- imap(models, ~ simulate_model_samples(.x, N, envs_template, Ci, PPFDs, param_means[[.y]]))
# result_arrays is a named list; each element is a 3D array [Ci, PPFD, N].
# Example access: result_arrays[["FvCB (1980)"]][ , "500", 1]

saveRDS(result_arrays, "result_arrays_scenario2.rds")