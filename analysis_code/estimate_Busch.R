library(tidyr)
library(dplyr)
library(ggpubr)

devtools::load_all()

# find Ko, kc, V_omax

calc_Anet <- function(Ci, KC, KO, Vomax){
  envs <- list(
    C_chl = Ci,
    O = rep(210, length(Ci)),
    PPFD = rep(1500, length(Ci)))

  pars <- list(
    R_d = 1.29,
    K_C = KC,
    K_O = KO,
    V_cmax = 104,
    V_omax = Vomax,
    J = 185,
    V_tpu = 11.5,
    max_alpha_G = 0.09,
    max_alpha_S = 0.38,
    N_max = 1.21
  )

  C_transition=NULL
  Ares <- Busch2018(envs=envs,pars=pars)
  return(Ares$An)
}


Ci <- c(50,75,150,225, 300, 387.5, 475, 562.5, 650, 750, 850, 1025, 1212.5, 1400, 1587.5, 1775)

An <- c(3, 5.5, 14.5, 22, 28, 32.5, 35, 35.5, 36.5, 37, 36.75, 36.5, 35.5, 35, 34.5, 34.25)

ACi_curve <- data.frame(Ci = Ci, A = An)

plot(ACi_curve$Ci, ACi_curve$A)

start_vals <- list(KO = 404, KC = 278, Vomax = 20)

fit <- nls(
  A ~ calc_Anet(Ci, KO, KC, Vomax),
  data = ACi_curve,
  start = start_vals,
  algorithm = "port",
  lower = c(20, 20, 1),
  upper = c(700,  800, 500),
  control = nls.control(maxiter = 500)
)


#############




library(pso)
objfunc <- function(x){
  Ci <- c(50,75,150,225, 300, 387.5, 475, 562.5, 650, 750, 850, 1025, 1212.5, 1400, 1587.5, 1775)
  envs <- list(
    C_chl = Ci,
    O = rep(210, length(Ci)),
    PPFD = rep(1500, length(Ci)))

  pars <- list(
    R_d = 1.29,
    K_C = x[1],
    K_O = x[2],
    V_cmax = 104,
    V_omax = x[3],
    J = 185,
    V_tpu = 11.5,
    max_alpha_G = 0.09,
    max_alpha_S = 0.38,
    N_max = 1.21
  )

  Asim <- Busch2018(envs=envs,pars=pars)$An
  Areal <- c(3, 5.5, 14.5, 22, 28, 32.5, 35, 35.5, 36.5, 37, 36.75, 36.5, 35.5, 35, 34.5, 34.25)

  fitness <- sqrt(mean((Areal - Asim)^2))
  return(fitness)
}


nvar <- 3
lb <- c(20, 20, 1)
ub <- c(700,  800, 500)

res <- psoptim(par = rep(NA,nvar), fn = objfunc,
               lower = lb, upper = ub,
               control = list(maxit = 5000, s = 100, maxit.stagnate = 500))

est.val <- res$par

simA <- calc_Anet(Ci, est.val[1], est.val[2], est.val[3])

plot(Ci, simA)

###############

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
