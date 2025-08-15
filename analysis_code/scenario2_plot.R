library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

devtools::load_all()

# Grid search for environmental input:
Ci <- c(20,50,100,150,200,300,400,500,600,700,800,1000,1200)
PPFDs <- c(20,100,300,500,700,1000,1200,1500)

envs_template <- list(
  C_i = Ci,
  O = rep(210, length(Ci))
)

# Define unique param_mean for each model
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

# Model registry: function + its parameter list
models <- list(
  `FvCB (1980)` = list(
    fun  = FvCB1980,
    pars = list(
      V_cmax=param_means[["FvCB (1980)"]]["Vcmax"], 
      J_max=param_means[["FvCB (1980)"]]["Jmax"], 
      R_d=param_means[["FvCB (1980)"]]["Rd"],
      K_C=param_means[["FvCB (1980)"]]["KC"], 
      K_O=param_means[["FvCB (1980)"]]["KO"], 
      gamma_star=param_means[["FvCB (1980)"]]["Gstar"],
      alpha_J=param_means[["FvCB (1980)"]]["alphaJ"],
      theta_J=param_means[["FvCB (1980)"]]["thetaJ"], 
      V_tpu=param_means[["FvCB (1980)"]]["Vtpu"]
    )
  ),
  `Harley (1992)` = list(
    fun  = Harley1992, 
    pars = list(
      V_cmax=param_means[["Harley (1992)"]]["Vcmax"], 
      J_max=param_means[["Harley (1992)"]]["Jmax"], 
      R_d=param_means[["Harley (1992)"]]["Rd"],
      K_C=param_means[["Harley (1992)"]]["KC"], 
      K_O=param_means[["Harley (1992)"]]["KO"], 
      gamma_star=param_means[["Harley (1992)"]]["Gstar"],
      gm=param_means[["Harley (1992)"]]["gm"]
    )
  ),
  `Ethier & Long (2004)` = list(
    fun  = EthierLong2004,
    pars = list(
      V_cmax=param_means[["Ethier & Long (2004)"]]["Vcmax"], 
      J_max=param_means[["Ethier & Long (2004)"]]["Jmax"], 
      R_d=param_means[["Ethier & Long (2004)"]]["Rd"],
      K_C=param_means[["Ethier & Long (2004)"]]["KC"], 
      K_O=param_means[["Ethier & Long (2004)"]]["KO"], 
      gamma_star=param_means[["Ethier & Long (2004)"]]["Gstar"],
      gm=param_means[["Ethier & Long (2004)"]]["gm"], 
      alpha_J=param_means[["Ethier & Long (2004)"]]["alphaJ"], 
      theta_J=param_means[["Ethier & Long (2004)"]]["thetaJ"]
    )
  ),
  `von Caemmerer (2000)` = list(
    fun  = Caemmerer2000,
    pars = list(
      V_cmax=param_means[["von Caemmerer (2000)"]]["Vcmax"], 
      J_max=param_means[["von Caemmerer (2000)"]]["Jmax"], 
      R_d=param_means[["von Caemmerer (2000)"]]["Rd"],
      K_C=param_means[["von Caemmerer (2000)"]]["KC"], 
      K_O=param_means[["von Caemmerer (2000)"]]["KO"], 
      gamma_star=param_means[["von Caemmerer (2000)"]]["Gstar"],
      alpha_J=param_means[["von Caemmerer (2000)"]]["alphaJ"], 
      theta_J=param_means[["von Caemmerer (2000)"]]["thetaJ"],
      V_tpu=param_means[["von Caemmerer (2000)"]]["Vtpu"], 
      alpha_tpu=param_means[["von Caemmerer (2000)"]]["atpu"]
    )
  ),
  `Yin (2004)` = list(
    fun  = Yin2004,
    pars = list(
      V_cmax=param_means[["Yin (2004)"]]["Vcmax"], 
      J_max=param_means[["Yin (2004)"]]["Jmax"], 
      R_d=param_means[["Yin (2004)"]]["Rd"],
      K_C=param_means[["Yin (2004)"]]["KC"], 
      K_O=param_means[["Yin (2004)"]]["KO"], 
      gamma_star=param_means[["Yin (2004)"]]["Gstar"],
      theta_J=param_means[["Yin (2004)"]]["thetaJ"], 
      phi2m=param_means[["Yin (2004)"]]["phi2"], 
      f_Q=param_means[["Yin (2004)"]]["fQ"],
      f_pseudo=param_means[["Yin (2004)"]]["fpseudo"], 
      f_cyc=param_means[["Yin (2004)"]]["fcyc"], 
      h=param_means[["Yin (2004)"]]["h"]
    )
  ),
  `Dubois (2007)` = list(
    fun  = Dubois2007,
    pars = list(
      V_cmax=param_means[["Dubois (2007)"]]["Vcmax"], 
      J_max=param_means[["Dubois (2007)"]]["J"], 
      R_d=param_means[["Dubois (2007)"]]["Rd"],
      K_C=param_means[["Dubois (2007)"]]["KC"], 
      K_O=param_means[["Dubois (2007)"]]["KO"], 
      gamma_star=param_means[["Dubois (2007)"]]["Gstar"]
    )
  )
)

run_model_once <- function(model_def, envs_base, PPFD_val){
  envs <- envs_base
  envs$PPFD <- rep(PPFD_val, length(envs$C_i))
  
  res <- model_def$fun(envs = envs, pars = model_def$pars)
  
  A <- res$An
  tibble(Ci = envs$C_i, A = A, PPFD = PPFD_val)
}

# Run all models × PPFDs:
df <- imap_dfr(models, function(mdef, mname){
  map_dfr(PPFDs, ~ run_model_once(mdef, envs_template, .x)) %>%
    mutate(Model = mname)
})

df$PPFD <- factor(df$PPFD, levels=PPFDs)

ggplot(df, aes(Ci, A, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~PPFD, ncol = 4) +
  labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
       y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
       color = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")