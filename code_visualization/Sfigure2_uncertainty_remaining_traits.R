library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(TeachingDemos)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")
load(file="results/Rdata/fig2_sampled_params.RData")

load(file = "results/Rdata/data_visualization/fig2_CIranges.RData")
bounds <- list(
  V_cmax     = c(50, 200),
  J_max      = c(100, 250),
  R_d        = c(0.5, 5),
  gamma_star = c(10, 35),
  V_tpu      = c(8, 20),
  gm         = c(0.1, 0.7),
  alpha_J    = c(0.3, 0.5),
  theta_J    = c(0.5, 0.9),
  K_C        = c(200, 300),
  K_O        = c(150, 300),
  alpha_G    = c(0.3, 0.8),
  phi2m      = c(0.5, 0.9),
  f_Q        = c(0.5, 1),
  f_pseudo   = c(0, 0.3),
  h          = c(3, 14/3),
  max_aG     = c(0.01, 0.2),
  max_aS     = c(0.01, 0.5),
  N_max      = c(0.1, 2),
  g_ch       = c(0.05, 0.7),
  g_wp       = c(0.1, 2),
  s          = c(0.3, 0.6),
  Phi2LL     = c(0.5, 1)
)

common_par <-c("theta_J","alpha_J", "K_C", "K_O") #,


newmodels <-  c("FvCB80", "Harley92","vonC.00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch17","Xiao21")


models <-  c("FvCB80", "Harley92","Caemmerer00",
             "Ethier04", "Yin04","Dubois07",
             "Tholen12", "Busch18","Xiao20")

df$data[df$data=="ACi(satQ)"] <- "Single A-Ci"
df$data[df$data=="ACi & AQ"] <- "A-Ci & A-Q"
df$data[df$data=="ACi(different Q)"] <- "Multiple A-Ci"
df <- df[df$data!="ACi(different Q) & CF",]

ind <- which(df$variable=="theta_J" & df$model=="Harley92")
df <- df[-ind,]

datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")
df$data <- factor(df$data, levels = datatypes)

df$model <- newmodels[match(df$model,models)]
df$model <- factor(df$model, levels = newmodels)

df <- df[df$variable%in%common_par,]
df$variable <- factor(df$variable, levels = common_par)

exp_vec <- c(   expression(theta[J]),
                expression(alpha[J]),
                expression(K[C]),
                expression(K[O]))


levels(df$variable)<- exp_vec

cols <- c("#1B9E77", "#D95F02", "#7570B3")
g1 <- ggplot(df,
       aes(x = data, y = range_perc, fill = data, colour= data)) +
  geom_boxplot( alpha = 0.5) +

  facet_grid(variable ~ model,
             switch = "y",
             scales = "free_y",
             labeller = labeller(variable = label_parsed)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 10),
    legend.position = "none",

    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "", y = "Relative CI width (%)",
       fill = "", colour= "") +

  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(0, 65)) +
  scale_y_continuous(breaks = seq(0, 60, 20))  # specify breaks
###################################################################################


load(file = "results/Rdata/data_visualization/fig2_CIranges.RData")

df$data[df$data=="ACi(satQ)"] <- "Single A-Ci"
df$data[df$data=="ACi & AQ"] <- "A-Ci & A-Q"
df$data[df$data=="ACi(different Q)"] <- "Multiple A-Ci"
df <- df[df$data!="ACi(different Q) & CF",]

ind <- which(df$variable=="theta_J" & df$model=="Harley92")
df <- df[-ind,]

datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")
df$data <- factor(df$data, levels = datatypes)

df$model <- newmodels[match(df$model,models)]
df$model <- factor(df$model, levels = newmodels)


vars <- c("V_cmax","J_max", "R_d", "gamma_star","theta_J", "alpha_J","K_C", "K_O",
          "gm","V_tpu","alpha_G", "phi2m", "f_Q","f_pseudo", "h","max_aG","max_aS","N_max",
          "g_ch",  "g_wp","s","Phi2LL" )

specific_par <- setdiff(vars,common_par)

df <- df[df$variable%in%specific_par,]
df$variable <- factor(df$variable, levels = specific_par)

df$model <- factor(df$model, levels = newmodels)


var_map <- c(
  V_cmax   = "V[cmax]",
  J_max    = "J[max]",
  R_d      = "R[d]",
  gamma_star = "Gamma^\"*\"",
  theta_J  = "theta[J]",
  alpha_J  = "alpha[J]",
  K_C      = "K[C]",
  K_O      = "K[O]",
  gm       = "g[m]",
  V_tpu    = "V[TPU]",
  alpha_G  = "alpha[G]",
  phi2m    = "Phi[IIm]",
  f_Q      = "f[Q]",
  f_pseudo = "f[pseudo]",
  h        = "h",
  g_ch     = "g[ch]",
  g_wp     = "g[wp]",
  max_aG   = "max[aG]",
  max_aS   = "max[aS]",
  N_max    = "N[max]",
  s        = "s",
  Phi2LL   = "Phi[IILL]"
)


combos <- tibble::tribble(
  ~model,     ~variable,
  "Harley92","gm",
  "vonC.00","V_tpu",
  "vonC.00","alpha_G",
  "Ethier04","gm",
  "Yin04","phi2m",
  "Yin04","f_Q",
  "Yin04","f_pseudo",
  "Yin04","h" ,
  "Dubois07","V_tpu",
  "Dubois07","alpha_G" ,
  "Tholen12","g_ch" ,
  "Tholen12","g_wp" ,
  "Busch17","V_tpu",
  "Busch17","max_aG",
  "Busch17","max_aS",
  "Busch17","N_max"  ,
  "Xiao21","gm",
  "Xiao21","s",
  "Xiao21","Phi2LL"
)

df2 <- semi_join(df, combos, by = c("model","variable")) %>%
  mutate(
    # make a model • variable label where variable is plotmath
    variable_lab = dplyr::recode(variable, !!!var_map),
    facet_lab = paste(model,"~", variable_lab),  # plain text + plotmath
    facet_lab = factor(
      facet_lab,
      levels = paste(combos$model,"~", dplyr::recode(combos$variable, !!!var_map))
    )
  )

df2$data <- factor(df2$data, levels = datatypes)
df2$model <- factor(df2$model, levels = newmodels)


cols <- c("#1B9E77", "#D95F02", "#7570B3")

g2 <- ggplot(df2, aes(y= range_perc, x = data, fill = data, colour = data)) +
  geom_boxplot( alpha = 0.50) +
  facet_wrap(~ facet_lab, scales = "free_y",
             labeller = label_parsed) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        strip.text.x = element_text(size=9),
        legend.text = element_text(size=11),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 1),
        panel.grid.major.x = element_blank()
  )+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols) +
  labs(x = "", y = "Relative CI width (%)",
       fill = "",  colour= "") +
  coord_cartesian(ylim = c(0, 65))+
  scale_y_continuous(breaks = seq(0, 60, 20))  # specify breaks

g2

library("ggpubr")
ggarrange(g1,g2,labels = c('a', 'b'), ncol=1, heights = c(0.8,1))

ggsave(filename = paste0("results/Figures/SFig2_remain.png.png"),width = 7, height = 8)


