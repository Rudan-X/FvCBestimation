library(tidyr)
library(dplyr)
library(ggpubr)

devtools::load_all()


calc_Anet <- function(Ci, Vcmax, J, Rd){
  envs <- list(
    C_chl = Ci,
    O = rep(210, length(Ci)))

  pars <- list(
    R_d= Rd,
    K_C= 400,
    K_O= 280,
    gamma_star= 45,
    V_cmax= Vcmax,
    J = J
  )

  C_transition=NULL
  Ares <- Dubois2007(envs=envs,pars=pars,Ct=C_transition)
  return(Ares$An)
}

calc_limit <- function(Ci, Vcmax, J, Rd){
  envs <- list(
    C_chl = Ci,
    O = rep(210, length(Ci)))

  pars <- list(
    R_d= Rd,
    K_C= 400,
    K_O= 280,
    gamma_star= 45,
    V_cmax= Vcmax,
    J = J
  )

  C_transition=NULL
  Ares <- Dubois2007(envs=envs,pars=pars,Ct=C_transition)
  return(Ares$Limitation)
}

C <- c(50,100,150,200,300,400,500,600,700,800,900,1000, 1200)
true.val <- c(150, 200, 1)
An <- calc_Anet(C, true.val[1], true.val[2], true.val[3])

ACi_curve <- data.frame(Ci = C, A = An)

plot(C, An)

start_vals <- list(Vcmax = 50, J = 100, Rd = 0.1)

fit <- nls(
  A ~ calc_Anet(Ci, Vcmax, J, Rd),
  data = ACi_curve,
  start = start_vals,
  algorithm = "port",
  lower = c(Vcmax = 0, J = 0, Rd = -3),
  upper = c(Vcmax = 500, J = 500, Rd = 50),
  control = nls.control(maxiter = 500)
)

est.val <- coef(fit)
true.Anet <- ACi_curve$An
sim.Anet <- calc_Anet(C, est.val[1], est.val[2], est.val[3])

true.Limit <- calc_limit(C, true.val[1], true.val[2], true.val[3])
sim.Limit <- calc_limit(C, est.val[1], est.val[2], est.val[3])


df <- data.frame(
  Ci = rep(C, 2),
  Anet = c(true.Anet, sim.Anet),
  Type = rep(c("True", "Simulated"), each = length(C)),
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
