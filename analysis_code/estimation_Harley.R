library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
devtools::load_all()

calc_true_Anet <- function(Ci, Vcmax, Jmax, Rd, gm){
  envs <- list(
    C_i = Ci,
    O = rep(210, length(Ci)),
    PPFD = rep(1500, length(Ci)))

  pars <- list(
    R_d = Rd,
    K_C = 274,
    K_O = 420,
    gamma_star = 43,
    V_cmax = Vcmax,
    J_max = Jmax
  )

  pars$gm <- gm

  Ares <- Harley1992(envs=envs,pars=pars,quadratic = TRUE,Ct = NULL)
  return(list(An=Ares$An, limit= Ares$Limitation))
}



sim_Anet <- function(A, Ci, Vcmax, Jmax, Rd){
  envs <- list(
    C_i = Ci,
    O = rep(210, length(Ci)),
    PPFD = rep(1500, length(Ci)))

  pars <- list(
    R_d = Rd,
    K_C = 274,
    K_O = 420,
    gamma_star = 43,
    V_cmax = Vcmax,
    J_max = Jmax
  )

  estim_gm <- estimate_gm_constantJ(A = A, Ci = Ci, pars = pars)

  envs$C_chl <- Ci - A/estim_gm
  Ares <- Harley1992(envs=envs,pars=pars,quadratic = FALSE,Ct = NULL)
  return(Ares$An)
}

sim_Limit <- function(A, Ci, Vcmax, Jmax, Rd){
  envs <- list(
    C_i = Ci,
    O = rep(210, length(Ci)),
    PPFD = rep(1500, length(Ci)))

  pars <- list(
    R_d = Rd,
    K_C = 274,
    K_O = 420,
    gamma_star = 43,
    V_cmax = Vcmax,
    J_max = Jmax
  )

  estim_gm <- estimate_gm_constantJ(A = A, Ci = Ci, pars = pars)

  envs$C_chl <- Ci - A/estim_gm
  Ares <- Harley1992(envs=envs,pars=pars,quadratic = FALSE,Ct = NULL)
  return(Ares$Limitation)
}


Ci <- c(20,50,100,150,200,300,400,500,600,700,800, 1000, 1200)

true.val <- c(70, 130, 0.5, 0.8)
true.sim <- calc_true_Anet(Ci, true.val[1], true.val[2], true.val[3], true.val[4])

plot(Ci, true.sim$An)

ACi_curve <- data.frame(Ci = Ci, A = true.sim$An)

start_vals <- list(Vcmax = 50, Jmax = 100, Rd = 0.2)

fit <- nls(
  A ~ sim_Anet(A, Ci, Vcmax, Jmax, Rd),
  data = ACi_curve,
  start = start_vals,
  algorithm = "port",
  lower = c(Vcmax = 0, Jmax = 0, Rd = -2),
  upper = c(Vcmax = 500, Jmax = 500, Rd = 50),
  control = nls.control(maxiter = 500)
)

est.val <- coef(fit)


true.Anet <- true.sim$An
true.Limit <- true.sim$limit

sim.Anet <- sim_Anet(ACi_curve$A, ACi_curve$Ci, est.val[1], est.val[2], est.val[3])

sim.Limit <- sim_Limit(ACi_curve$A, ACi_curve$Ci, est.val[1], est.val[2], est.val[3])


df <- data.frame(
  Ci = rep(Ci, 2),
  Anet = c(true.Anet, sim.Anet),
  Type = rep(c("True", "Simulated"), each = length(Ci)),
  Limitation = c(true.Limit, sim.Limit)
)

ggplot(df, aes(x = Ci, y = Anet, color = Type)) +
  geom_line(aes(linetype = Type), size = 1) +
  geom_point(aes(shape = Limitation), size = 3, alpha=0.7) +
  scale_linetype_manual(values = c("True" = "dashed", "Simulated" = "dashed")) +
  scale_color_manual(values = c("True" = "blue", "Simulated" = "red")) +
  scale_shape_manual(values = c("Wc" = 15, "Wj" = 16, "Wp" = 17)) +  # You can change shapes if needed
  labs(
    x = expression(C[i]~"(µmol mol"^{-1}*")"),
    y = expression(A[n]~"(µmol m"^{-2}~s^{-1}*")"),
    color = "Curve",
    linetype = "Curve",
    shape = "Limitation"
  ) +
  theme_minimal(base_size = 14)

