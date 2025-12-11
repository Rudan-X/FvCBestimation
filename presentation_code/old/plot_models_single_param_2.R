library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

devtools::load_all()
# check check
# Grid search for environmental input:
# Ci <- c(20,50,100,150,200,300,400,500,600,700,800,1000,1200)
# PPFDs <- c(20,100,300,500,700,1000,1200,1500)
Ci     <- c(50,100,150,200,300,500,600,700,800,1000,1200)
PPFDs  <- c(50,150,300,500,1100,1500,1800)

envs_template <- list(
  C_i = Ci,
  O = rep(210, length(Ci))
)

Dubois_param <- c(Vcmax = 166, Jmax= 239, Rd=4.5)
param <- c(
  Vcmax = 104,
  Jmax = 200,
  Rd = 1.29,
  KC = 270,
  KO = 165,
  Gstar = 34.2,
  alphaJ = 0.3,
  thetaJ = 0.8,
  Vtpu = 11.5,
  gm = 0.3,
  atpu = 0.3,
  phi2 = 0.85,
  fQ = 1,
  fpseudo = 0,
  fcyc = NULL,
  h = 14/3,
  gwp = 1.35,
  gch = 0.34,
  Sco = 3070,
  max_ag = 0.09,
  max_as = 0.38,
  Nmax = 1.21
)
# Model registry: function + its parameter list
models <- list(
  list(
    fun  = FvCB1980,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      alpha_J=param["alphaJ"],
      theta_J=param["thetaJ"], V_tpu=param["Vtpu"]
    )
  ),
  list(
    fun  = Harley1992,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      gm=param["gm"]
    )
  ),
  list(
    fun  = Caemmerer2000,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      alpha_J=param["alphaJ"], theta_J=param["thetaJ"],
      V_tpu=param["Vtpu"],alpha_tpu=param["atpu"]
    )
  ),
  list(
    fun  = EthierLong2004,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      gm=param["gm"], alpha_J=param["alphaJ"], theta_J=param["thetaJ"]
    )
  ),
  list(
    fun  = Yin2004,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      theta_J=param["thetaJ"], phi2m = param["phi2"],f_Q = param["fQ"],
      f_pseudo = param["fpseudo"],f_cyc = param["fcyc"],h = param["h"]
    )
  ),
  list(
    fun  = Dubois2007,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      alpha_J=param["alphaJ"],theta_J=param["thetaJ"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      V_tpu=param["Vtpu"],alpha_tpu=param["atpu"]
    )
  ),
  list(
    fun  = Tholen2012,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      K_C=param["KC"], K_O=param["KO"],  gamma_star=param["Gstar"],
      g_wp=param["gwp"], g_ch=param["gch"], alpha_J=param["alphaJ"], theta_J=param["thetaJ"]
    )
  ),
  list(
    fun  = Busch2018,
    pars = list(
      V_cmax=param["Vcmax"], J_max=param["Jmax"], R_d=param["Rd"],
      alpha_J=param["alphaJ"],theta_J=param["thetaJ"],
      K_C=param["KC"], K_O=param["KO"], gamma_star=param["Gstar"],
      V_tpu=param["Vtpu"],max_alpha_G=param["max_ag"], max_alpha_S=param["max_as"],
      N_max=param["Nmax"]
    )
  )
)



model_names <- c("FvCB (1980) (Wc,Wj) (Ci~=Cc)", "Harley (1992) (Ac,Aj) (Ci=Cc-A/gm)",
                 "von Caemmerer (2000) (Wc,Wj,Wp) (Ci~=Cc)","Ethier & Long (2004) (Ac,Aj) (Ci=Cc-A/gm)",
                 "Yin (2004) (Wc,Wj) (Ci~=Cc)","Dubois (2007) (Ac,Aj,Ap) (Ci~=Cc)",
                 "Tholen (2012) (Wc,Wj) (Ci=Cc-A/gm)","Busch (2018) (Wc,Wj,Wp) (Ci~=Cc)")

names(models) <- model_names
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

new_names <- c("FvCB (1980) (Wc,Wj) (Ci~=Cc)", "Yin (2004) (Wc,Wj) (Ci~=Cc)",
               "Harley (1992) (Ac,Aj) (Ci=Cc-A/gm)", "Ethier & Long (2004) (Ac,Aj) (Ci=Cc-A/gm)",
               "Tholen (2012) (Wc,Wj) (Ci=Cc-A/gm)",
               "von Caemmerer (2000) (Wc,Wj,Wp) (Ci~=Cc)","Dubois (2007) (Ac,Aj,Ap) (Ci~=Cc)",
               "Busch (2018) (Wc,Wj,Wp) (Ci~=Cc)")

df$PPFD <- factor(df$PPFD, levels=PPFDs)
df$Model <- factor(df$Model, levels= new_names)
ggplot(df, aes(Ci, A, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~PPFD, ncol = 4, dir = "v") +
  labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
       y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
       color = "PPFD") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ggsave("results/Figures/Single_parameter.png",width = 12, height = 7)
