library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
devtools::load_all()

Ci <- c(20,50,100,150,200,300,400,500,600,700,800, 1000, 1200)
envs <- list(
  C_i=Ci,
  PPFD=rep(1800,length(Ci)),
  O=rep(210,length(Ci))) # , )

# Base parameter list
base_pars <- list(
  R_d = 0.5,
  K_C = 274,
  K_O = 420,
  gm = 0.8,
  gamma_star = 43,
  V_cmax = 70,
  J_max = 130
)

# gm values to test
gm_values <- c(0.8, 0.4, 0.2, 0.1, 0.05)

# Run model for each gm
df_all <- lapply(gm_values, function(gm_val) {
  pars <- base_pars
  pars$gm <- gm_val
  res <- Harley1992(envs, pars, quadratic = TRUE)
  data.frame(Ci = envs$C_i, Anet = res$An, gm = gm_val)
}) %>% bind_rows()

# Plot with ggplot
ggplot(df_all, aes(x = Ci, y = Anet, color = factor(gm))) +
  geom_line(size = 1) +
  geom_point() +
  labs(
    x = expression(C[i]~"(µmol mol"^-1*")"),
    y = expression(A[net]~"(µmol m"^-2~s^-1*")"),
    color = expression(g[m]~"(mol m"^-2~s^-1~bar^-1*")")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Anet=Harley1992(envs,pars)
#
# df <- as.data.frame(Anet)
# df$Ci <- envs$C_i # add a time or observation index
# df <- df[,-4]
#
# # Reshape to long format
# df_long <- df %>%
#   pivot_longer(cols = -Ci, names_to = "Limitation", values_to = "Value")
#
#
# df_long <- df_long %>%
#   mutate(Linetype = ifelse(Limitation == "An", "solid", "dashed"))
#
# df_long$Limitation <- factor(df_long$Limitation, levels=c("Ac","Aj","An"))
#
#
# ggplot(df_long, aes(x = Ci, y = Value, color = Limitation, linetype = Linetype)) +
#   geom_line(size = 1, alpha = 0.9) +
#   scale_linetype_identity() +  # Use values from the Linetype column directly
#   theme_minimal(base_size = 14) +
#   labs(
#     x = expression(C[i]~"(µmol mol"^{-1}*")"),
#     y = expression(An~"(µmol m"^{-2}~s^{-1}*")"),
#     color = NULL
#   ) +
#   theme(
#     legend.position = "right"
#   )
