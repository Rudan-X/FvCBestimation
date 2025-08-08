

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

devtools::load_all()


Ci <- c(20,50,100,150,200,300,400,500,600,700,800, 1000, 1200)
envs <- list(
  C_i = Ci,
  PPFD = rep(1500,length(Ci)),
  O = rep(210,length(Ci)))

pars <- list(
  R_d=2,
  K_C=268,
  K_O=165,
  gamma_star=38,
  V_cmax=70,
  J_max=130,
  alpha_J=0.33,
  theta_J=0.825,
  V_tpu=9
)

A_1980=FvCB1980(envs,pars)


pars <- list(
  V_cmax = 70,
  J_max = 130,
  K_C = 274,
  K_O = 420,
  R_d = 0.5,
  gamma_star = 43,
  gm = 0.5
)
A_1992 <- Harley1992(envs=envs,pars=pars)


pars <- list(
  V_cmax = 70,
  J_max = 130,
  K_C = 274,
  K_O = 420,
  R_d = 0.5,
  gamma_star = 43,
  gm = 0.5,
  alpha_J=0.33,
  theta_J=0.825
)

A_2004 <- EthierLong2004(envs=envs,pars=pars)


pars <- list(
  V_cmax = 70,
  J_max = 130,
  K_C = 274,
  K_O = 420,
  R_d = 0.5,
  gamma_star = 43,
  alpha_J=0.33,
  theta_J=0.825,
  V_tpu = 9,
  alpha_tpu= 0.3
)

A_2000 <- Caemmerer2000(envs=envs,pars=pars)


pars <- list(
  V_cmax = 70,
  J = 130,
  K_C = 274,
  K_O = 420,
  R_d = 0.5,
  gamma_star = 43
)

A_2007 <- Dubois2007(envs=envs,pars=pars)


pars <- list(
  R_d = 1.29,
  K_C = 404,
  K_O = 278,
  gamma_star = 43,
  V_cmax = 104,
  J_max = 185,
  theta_J = 0.7,
  phi2m = 0.85,
  f_Q = 1,
  f_pseudo = 0,
  f_cyc = NULL,
  h = 14
)

A_2004 <- Yin2004(envs=envs,pars=pars)
#
# PARs <- c(10,20,50,70,100,300,500,700,1000,1200, 1500)
#
# Ccs <- c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200)
# ACi_PAR <- matrix(nrow = length(PARs),ncol  = length(Ccs))
# rownames(ACi_PAR) <- as.character(PARs)
# colnames(ACi_PAR) <- as.character(Ccs)
#
# ACi_PAR_lim <- ACi_PAR
# for (i in 1:length(PARs)){
#   PAR <- rep(PARs[i],length(Ccs))
#   env <- list(C_chl= Ccs, PPFD=PAR)
#   parameters <- list(
#     V_cmax=x[1],
#     J_max=x[2],
#     V_tpu=x[3],
#     K_CO=x[4],
#     R_d=x[5],
#     gamma_star=x[6],
#     phi_J=0.33,
#     theta_J=0.825
#   )
#   res <- FvCB1980(envs = env, pars = parameters)
#   ACi_PAR[i,] <- res$An
#   ACi_PAR_lim[i,] <- res$Limitation
# }
