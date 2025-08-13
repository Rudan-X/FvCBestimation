# nolint start: object_name_linter
library(dplyr)
library(purrr)
library(devtools)
# Set working directory
setwd("C:\\Users\\Salma\\Documents\\Salma_Project_Data_All\\Kaib_paper\\C3_primary_model_julia_All_Mod\\FvCBestimation-master")


# Load your package functions and models
devtools::load_all()
models <- load_FvCBmodels()

# Grids for environment
Ci    <- c(20,50,100,150,200,300,400,500,600,700,800,1000,1200)
PPFDs <- c(20,100,300,500,700,1000,1200,1500)

envs_template <- list(
  C_i = Ci,
  O   = rep(210, length(Ci))
)

# Base parameters (means)
param_mean <- c(
  Vcmax   = 70,
  Jmax    = 130,
  J       = 130,
  Rd      = 2,
  KC      = 270,
  KO      = 165,
  Gstar   = 38,
  alphaJ  = 0.3,
  thetaJ  = 0.8,
  Vtpu    = 9,
  gm      = 0.3,
  atpu    = 0.3,
  phi2m   = 0.85,
  fQ      = 1,
  fpseudo = 0,
  fcyc    = NA_real_,  # set later if needed
  h       = 14
)

# Function to generate parameter samples for all parameters at once (samples x params)
sample_params_matrix <- function(param_mean, pct_range = 0.1, N) {
  samples <- sapply(param_mean, function(mu) {
    lo <- mu * (1 - pct_range)
    hi <- mu * (1 + pct_range)
    runif(N, min = lo, max = hi)
  })
  
  # Ensure 'fcyc' is NA_real_ for all samples (if needed)
  samples[, "fcyc"] <- NA_real_
  
  # Return matrix: rows = samples, cols = parameters
  return(samples)
}

# Function to simulate all models using the SAME parameter samples
simulate_all_models_same_samples <- function(models, param_samples, envs_template, Ci, PPFDs) {
  N <- nrow(param_samples)
  
  result_arrays <- list()
  
  for (model_name in names(models)) {
    model_def <- models[[model_name]]
    n_ci <- length(Ci)
    n_ppfd <- length(PPFDs)
    
    # Initialize output array: Ci x PPFD x samples
    out <- array(NA_real_, dim = c(n_ci, n_ppfd, N),
                 dimnames = list(
                   Ci = as.character(Ci),
                   PPFD = as.character(PPFDs),
                   sample = paste0("s", seq_len(N))
                 ))
    
    for (s in seq_len(N)) {
      # Extract parameter vector for sample s
      p <- param_samples[s, ]
      
      # Map parameters for the model
      pars <- model_def$map(p)
      
      for (j in seq_along(PPFDs)) {
        envs <- envs_template
        envs$PPFD <- rep(PPFDs[j], length(Ci))
        
        res <- model_def$fun(envs = envs, pars = pars)
        A <- res$An
        
        if (length(A) != n_ci) {
          stop("Model returned An length ", length(A), " but expected ", n_ci)
        }
        
        out[, j, s] <- A
      }
    }
    
    result_arrays[[model_name]] <- out
  }
  
  return(result_arrays)
}

# ------------------ Run scenarios for different parameter ranges ------------------

N <- 10000                # number of parameter samples
pct_ranges <- c(0.10, 0.20, 0.30, 0.50)  # parameter uncertainty ranges

scenario_samples <- list()
scenario_results <- list()

for (pct in pct_ranges) {
  cat("Sampling parameters with ±", pct*100, "% uncertainty...\n")
  samples <- sample_params_matrix(param_mean, pct_range = pct, N = N)
  
  cat("Simulating all models with these samples...\n")
  sim_results <- simulate_all_models_same_samples(models, samples, envs_template, Ci, PPFDs)
  
  # Store with descriptive names
  samples_name <- paste0("scenario_", pct*100, "perc_samples")
  results_name <- paste0("scenario_", pct*100, "perc_results")
  
  scenario_samples[[samples_name]] <- samples
  scenario_results[[results_name]] <- sim_results
}

# Save results if desired:
saveRDS(scenario_samples, "scenario_samples.rds")
saveRDS(scenario_results, "scenario_results.rds")
