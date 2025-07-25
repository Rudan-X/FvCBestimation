library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(pso)
library(devtools)

devtools::load_all()
calc_Anet <- function(x){

  PARs <- c(10,20,50,70,100,300,500,700,1000,1200)

  Ccs <- c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200)
  ACi_PAR <- matrix(nrow = length(PARs),ncol  = length(Ccs))
  rownames(ACi_PAR) <- as.character(PARs)
  colnames(ACi_PAR) <- as.character(Ccs)

  ACi_PAR_lim <- ACi_PAR
  for (i in 1:length(PARs)){
    PAR <- rep(PARs[i],length(Ccs))
    env <- list(C_chl= Ccs, PPFD=PAR)
    parameters <- list(
      V_cmax=x[1],
      J_max=x[2],
      V_tpu=x[3],
      K_CO=x[4],
      R_d=x[5],
      gamma_star=x[6],
      phi_J=0.33,
      theta_J=0.825
    )
    res <- FvCB1980(envs = env, pars = parameters, Ct = NULL)
    ACi_PAR[i,] <- res$An
    ACi_PAR_lim[i,] <- res$Limitation
  }

  return(list(An=ACi_PAR, Limitation=ACi_PAR_lim))
}

true.val <- c(70, 130,9,700, 2, 38)
ACi <- list(A= calc_Anet(true.val)$An)

#### parameter estimation ####
# Objective function
objfunc <- function(x){
  Anet <- calc_Anet(x)$An
  fitness <- sqrt(mean((ACi$A - Anet)^2))
  return(fitness)
}


nvar <- 6
lb <- seq(0.01, nvar)
ub <- c(300, 500, 20, 2000, 20, 100)

res <- psoptim(par = rep(NA,nvar), fn = objfunc,
               lower = lb, upper = ub,
               control = list(maxit = 5000, s = 100, maxit.stagnate = 500))



true.Anet=calc_Anet(true.val)$An#
sim.Anet=calc_Anet(res$par)$An

true.Limit=calc_Anet(true.val)$Limitation#
sim.Limit=calc_Anet(res$par)$Limitation

df <- melt(true.Anet)
colnames(df) <- c("PAR","Ci","Anet")
df$PAR <- as.factor(df$PAR)

# Plot
ggplot(df, aes(x = Ci, y = Anet, color = PAR)) +
  geom_point() +
  geom_line() +
  labs(x = expression(Ci~(mu*mol/mol)),
       y = expression(An~(mu*mol~m^-2~s^-1)),
       color = "PAR level") +
  theme_minimal() +
  theme(legend.position = "right")





#############################
library(tidyverse)

# Assume true.Anet and sim.Anet are matrices
# Rows = PAR levels, Cols = Ci levels
PAR_vals <- c(10,20,50,70,100,300,500,700,1000,1200)
Ci_vals <- c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200)


# Convert to tidy data frame
df_true <- as.data.frame(true.Anet) %>%
  mutate(PAR = PAR_vals, Type = "True") %>%
  pivot_longer(-c(PAR, Type), names_to = "Ci", values_to = "Anet") %>%
  mutate(Ci = as.numeric(Ci))

df_sim <- as.data.frame(sim.Anet) %>%
  mutate(PAR = PAR_vals, Type = "Simulated") %>%
  pivot_longer(-c(PAR, Type), names_to = "Ci", values_to = "Anet") %>%
  mutate(Ci = as.numeric(Ci))

# Combine both
df_all <- bind_rows(df_true, df_sim)

ggplot(df_all, aes(x = Ci, y = Anet, color = factor(PAR),
                   shape = Type, linetype = Type)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = expression(C[i]~(mu*mol~mol^{-1})),
    y = expression(A[n]~(mu*mol~m^{-2}~s^{-1})),
    color = "PAR level",
    shape = "Data type",
    linetype = "Data type"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
#############################

true.val - res$par
