library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

devtools::load_all()

Ci <- c(20,50,100,150,200,300,400,500,600,700,800, 1000, 1200, 1500  )

envs <- list(
  C_i = Ci,
  O = rep(210, length(Ci)),
  PPFD = rep(1500, length(Ci)))

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


Anet=Yin2004(envs,pars)

df <- as.data.frame(Anet)
df$Ci <- envs$C_i # add a time or observation index
df <- df[,c("Ac","Aj","An","Ci")]

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = -Ci, names_to = "Limitation", values_to = "Value")


df_long <- df_long %>%
  mutate(Linetype = ifelse(Limitation == "An", "solid", "dashed"))

df_long$Limitation <- factor(df_long$Limitation, levels=c("Ac","Aj","Ap","An"))
# Plot
ggplot(df_long, aes(x = Ci, y = Value, color = Limitation, linetype = Linetype)) +
  geom_line(size = 1, alpha = 0.9) +
  scale_linetype_identity() +  # Use values from the Linetype column directly
  theme_minimal(base_size = 14) +
  labs(
    x = expression(C[i]~"(µmol mol"^{-1}*")"),
    y = expression(An~"(µmol m"^{-2}~s^{-1}*")"),
    color = NULL
  ) +
  theme(
    legend.position = "right"
  )
