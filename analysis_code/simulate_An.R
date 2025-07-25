

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

devtools::load_all()

envs <- list(
  C_chl=c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200),
  PPFD=rep(1800,14)) # ,  O=0.21)

pars <- list(
  R_d=2,
  K_CO=700,
  # K_C=268,
  # K_O=165, #084,
  gamma_star=38,
  V_cmax=70,
  J_max=130,
  phi_J=0.33,
  theta_J=0.825,
  V_tpu=9
)


gs <- simulate_FvCB1980(envs,pars)

ggarrange(gs$Wmin, gs$An, ncol = 1, labels = c("a", "b"))
