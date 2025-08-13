library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

devtools::load_all()
# check check
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
            h = 14,
            gwp = 1.35,
            gch = 0.34,
            Sco = 2592,
            max_ag = 0.09,
            max_as = 0.38,
            Nmax = 1.21
           )
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
      theta_J=param["thetaJ"], phi2m = param["phi2"],f_Q = param["fQ"],
      f_pseudo = param["fpseudo"],f_cyc = param["fcyc"],h = param["h"]
    )
  ),
  `Dubois (2007)` = list(
    fun  = Dubois2007,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      alpha_J=param["alphaJ"],theta_J=param["thetaJ"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      V_tpu=param["Vtpu"],alpha_tpu=param["atpu"]
    )
  ),
  `Tholen (2012)` = list(
    fun  = Tholen2012,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], S_co=param["Sco"],
      g_wp=param["gwp"], g_ch=param["gch"], alpha_J=param["alphaJ"], theta_J=param["thetaJ"]
    )
  ),
  `Busch (2018)` = list(
    fun  = Busch2018,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      alpha_J=param["alphaJ"],theta_J=param["thetaJ"],
      K_C=param["KC"], K_O=param["KO"], S_co=param["Sco"],
      V_tpu=param["Vtpu"],max_alpha_G=param["max_ag"], max_alpha_S=param["max_as"],
      N_max=param["Nmax"]
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
df$Model <- factor(df$Model, levels= c("FvCB (1980)", "Harley (1992)","von Caemmerer (2000)","Ethier & Long (2004)",
                                       "Yin (2004)","Dubois (2007)","Tholen (2012)","Busch (2018)"))
ggplot(df, aes(Ci, A, color = PPFD)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~Model, ncol = 4) +
  labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
       y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
       color = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

check=data.frame(dA=df$A[df$Model=="FvCB (1980)"]-df$A[df$Model=="Dubois (2007)"],
                 Ci=df$Ci[df$Model=="FvCB (1980)"], PPFD=df$PPFD[df$Model=="FvCB (1980)"])
# df <- read.csv("../FvCBestimation_not_uploaded/data/2022_ACi_rawData_barley.csv")
#
# df <- df[21:30,1:2]
#
# plot(df$Ci, df$Photo)
