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

param <- c(Vcmax = 70,
                Jmax = 130,
                J = 130,
                Rd = 2,
                KC = 270,
                KO = 165,
                Gstar = 38,
                alphaJ = 0.3,
                thetaJ = 0.8,
                Vtpu = 9,
                gm = 0.3,
                atpu = 0.3,
                phi2 = 0.85,
                fQ = 1,
                fpseudo = 0,
                fcyc = NULL,
                h = 14)
# Model registry: function + its parameter list
models <- list(
  `FvCB (1980)` = list(
    fun  = FvCB1980,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      alpha_J=param["alphaJ"],
      theta_J=param["thetaJ"], V_tpu=param["Vtpu"]
    )
  ),
  `Harley (1992)` = list(
    fun  = Harley1992,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      gm=param["gm"]
    )
  ),
  `Ethier & Long (2004)` = list(
    fun  = EthierLong2004,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      gm=param["gm"], alpha_J=param["alphaJ"], theta_J=param["thetaJ"]
    )
  ),

  `von Caemmerer (2000)` = list(
    fun  = Caemmerer2000,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      alpha_J=param["alphaJ"], theta_J=param["thetaJ"],
      V_tpu=param["Vtpu"],alpha_tpu=param["atpu"]
    )
  ),
  `Yin (2004)` = list(
    fun  = Yin2004,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      theta_J=param["thetaJ"], phi2m = param["phi2m"],f_Q = param["fQ"],
      f_pseudo = param["fpseudo"],f_cyc = param["fcyc"],h = param["h"]
    )
  ),
  `Dubois (2007)` = list(
    fun  = Dubois2007,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["J"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"]
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
