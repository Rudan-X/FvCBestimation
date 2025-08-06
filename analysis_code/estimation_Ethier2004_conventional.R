# Load required package
library(ggplot2)
devtools::load_all()
# Constants from Harley et al. 1992 (25 °C, 21% O2)

##############################
# simulate ideal curve:
Cc <- c(50,100,150,200,300,400,500,600,700,800, 1000, 1200)
Oc <- 210
K_C <- 274
K_O <-  420
envs <- list(
  C_chl = Cc,
  O = rep(Oc, length(Cc)),
  PPFD = rep(1500, length(Cc)))

pars <- list(
  R_d= 0.5,
  K_CO=K_C*(1+Oc/K_O),
  gamma_star= 45,
  V_cmax= 56,
  J_max = 112,
  phi_J= 0.18,
  theta_J= 0.9
)

C_transition=NULL
Ares <- Ethier2004(envs=envs,pars=pars,Ct=C_transition)

ideal_curve <- data.frame(A = Ares$An, Ac= Ares$Ac, Aj= Ares$Aj, Cc = Cc)

gm <- 0.2
ideal_curve$Ci <- ideal_curve$Cc + ideal_curve$A/gm

plot(ideal_curve$Ci, ideal_curve$A)
##############################
calc_Anet <- function(Ci, Vcmax, Jmax, gamma_s, KmCO){ #, Rd, KO
  envs <- list(
    C_chl = Ci,
    O = rep(210, length(Ci)),
    PPFD = rep(1500, length(Cc)))

  pars <- list(
    R_d= 0.5,
    K_CO=KmCO,
    gamma_star= gamma_s,
    V_cmax= Vcmax,
    J_max = Jmax,
    phi_J= 0.18,
    theta_J= 0.9
  )

  C_transition=NULL
  Ares <- Ethier2004(envs=envs,pars=pars,Ct=C_transition)
  return(Ares$An)
}

start_vals <- list(Vcmax = 20, Jmax = 50, gamma_s = 20, KmCO = 200)

# Fit model
fit <- nls(
  A ~ calc_Anet(Ci, Vcmax, Jmax, gamma_s, KmCO),
  data = ideal_curve,
  start = start_vals,
  algorithm = "port",
  lower = c(Vcmax = 5, Jmax = 5, gamma_s = 5, KmCO = 5),
  upper = c(Vcmax = 500,Jmax = 500, gamma_s = 100, KmCO = 700),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

estim <- coef(fit)
summary(fit)


envs <- list(
  C_chl = ideal_curve$Ci ,
  O = rep(Oc, length(Cc)),
  PPFD = rep(1500, length(Cc)))

pars <- list(
  R_d= 0.5,
  K_CO=estim["KmCO"],
  gamma_star= estim["gamma_s"],
  V_cmax= estim["Vcmax"],
  J_max = estim["Jmax"],
  phi_J= 0.18,
  theta_J= 0.9
)

C_transition=NULL
Ares <- Ethier2004(envs=envs,pars=pars,Ct=C_transition)

fitted_curve <- data.frame(A = Ares$An, Ac= Ares$Ac, Aj=Ares$Aj, Cc = Cc)

gm <- 0.2

fitted_curve$Ci <- fitted_curve$Cc + fitted_curve$A/gm
plot(fitted_curve$Ci, fitted_curve$Aj)
