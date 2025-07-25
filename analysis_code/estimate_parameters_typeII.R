library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(pso)

devtools::load_all()
calc_Anet <- function(x){
  parameters <- list(
    V_cmax=x[1],
    J_max=x[2],
    V_tpu=x[3],
    K_CO=x[4],
    R_d=2, # x[5],
    gamma_star=38, #x[6],
    phi_J=0.33,
    theta_J=0.825
  )
  env.setting <- list(
    C_chl=c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200),
    PPFD=rep(1800,14))

  C_transition=NULL
  Ares <- FvCB1980(envs=env.setting,parameters,Ct=C_transition)

  return(list(An=Ares$An, Limitation=Ares$Limitation))
}

true.val <- c(70, 130,9,700)
ACi <- list(A= calc_Anet(true.val)$An)

#### parameter estimation ####
# Objective function
objfunc <- function(x){

  Anet <- calc_Anet(x)$An
  fitness <- sqrt(mean((ACi$A - Anet)^2))
  return(fitness)
}


nvar <- 4
lb <- seq(0.01, nvar)
ub <- c(300,500,20, 2000)

res <- psoptim(par = rep(NA,nvar), fn = objfunc,
               lower = lb, upper = ub,
               control = list(maxit = 5000, s = 100, maxit.stagnate = 500))



true.Anet=calc_Anet(true.val)$An#
sim.Anet=calc_Anet(res$par)$An

true.Limit=calc_Anet(true.val)$Limitation#
sim.Limit=calc_Anet(res$par)$Limitation

envs <- list(
  C_chl=c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200),
  PPFD=1800)

Ci_vals   <- envs$C_chl

df <- data.frame(
  Ci = rep(Ci_vals, 2),
  Anet = c(true.Anet, sim.Anet),
  Type = rep(c("True", "Simulated"), each = length(Ci_vals)),
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


true.val-res$par
