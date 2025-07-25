library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(pso)
library(reshape2)
library(broom)

devtools::load_all()

# low PARs from 10 to 100 are essential to capture the more points in slope-intercept regression
# This aligns with studies Laisk 1977; Brooks & Farquhar 1985 saying about sub-saturating irradiances
PARs <- c(10,20,50,70,100,300,500,700)

Ccs <- c(20,50,100,150,200,300,400,500,600,700,800,900,1000, 1200)
ACi_PAR <- matrix(nrow = length(PARs),ncol  = length(Ccs))
for (i in 1:length(PARs)){
  PAR <- rep(PARs[i],length(Ccs))
  envs <- list(
    C_chl= Ccs,
    PPFD=PAR) # ,  O=0.21)
  pars <- list(
    V_cmax=70,
    J_max=130,
    V_tpu=9,
    K_CO=700,
    R_d=2,
    gamma_star=38,
    phi_J=0.33,
    theta_J=0.825
  )
  ACi_PAR[i,] <- FvCB1980(envs,pars)$An
}


rownames(ACi_PAR) <- as.character(PARs)
colnames(ACi_PAR) <- as.character(Ccs)

df <- melt(ACi_PAR)
colnames(df) <- c("PAR","Ci","Anet")
df$PAR <- as.factor(df$PAR)

# Plot
ggplot(df, aes(x = Ci, y = Anet, color = PAR, linetype = PAR)) +
  geom_point() +
  geom_line() +
  labs(x = expression(Ci~(mu*mol/mol)),
       y = expression(An~(mu*mol~m^-2~s^-1)),
       color = "PAR level") +
  theme_minimal() +
  theme(legend.position = "right")



# Slope-intercept regression to find Rd and gammma*
# doi: 10.1111/pce.12562
low_Ci_threshold <- 100

# Fit linear models per PAR level
fits <- df %>%
  filter(Ci <= low_Ci_threshold) %>%
  group_by(PAR) %>%
  do(tidy(lm(Anet ~ Ci, data = .))) %>%
  select(PAR, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(slope = Ci, intercept = `(Intercept)`)

# Step 2: Slope-intercept regression: b ~ m
slope_intercept_fit <- lm(intercept ~ slope, data = fits)
summary(slope_intercept_fit)

# Step 3: Extract Ci* and Rd
m_reg <- coef(slope_intercept_fit)[["slope"]]
b_reg <- coef(slope_intercept_fit)[["(Intercept)"]]

Ci_star <- -m_reg  # Ci* ~= gammma* == -slope
Rd <- b_reg        # Rd == intercept

cat("Estimated Ci* =", Ci_star, "\n")
cat("Estimated Rd =", Rd, "\n")

# Optional: plot slope-intercept regression
ggplot(fits, aes(x = slope, y = intercept)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Slope (m)", y = "Intercept (b)",
       title = "Slope–Intercept Regression (b ~ m)") +
  annotate("text", x = min(fits$slope), y = max(fits$intercept),
           label = paste0("Ci* = ", round(Ci_star, 2),
                          "\nRd = ", round(Rd, 2)),
           hjust = 0, vjust = 1)
