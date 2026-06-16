library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(ggpubr)

load(file="results/Rdata/fig2_sampled_params.RData")

true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")[,-1]

true_df <- melt(true_params_df)

true_df$Models <- "All models"
true_df$Models[true_df$variable%in%c("gm")] <- "Harley92,Ethier04&Xiao21"# "Harley92,Ethier04&Xiao20"
true_df$Models[true_df$variable%in%c("V_tpu")] <- "vonC00,Dubois07&Busch18"# "vonCaemmerer00,Dubois07&Busch18"
true_df$Models[true_df$variable%in%c("alpha_G")] <- "vonC00&Dubois07"# "vonCaemmerer00,Dubois07"
true_df$Models[true_df$variable%in%c("g_ch","g_wp")] <- "Tholen12" # "Tholen12"

true_df$Models[true_df$variable%in%c("h","f_Q","f_pseudo","phi2m")] <- "Yin04"# "Yin04"
true_df$Models[true_df$variable%in%c("max_aG","max_aS", "N_max")] <- "Busch18"#"Busch18"
true_df$Models[true_df$variable%in%c("s","Phi2LL")] <- "Xiao21"# "Xiao20"

vars <- c("V_cmax","J_max", "R_d", "gamma_star","theta_J", "alpha_J","K_C", "K_O",
          "gm","V_tpu","alpha_G", "phi2m", "f_Q","f_pseudo", "h","max_aG","max_aS","N_max",
          "g_ch",  "g_wp","s","Phi2LL" )

true_df$variable <- factor(true_df$variable, levels=vars)

true_df$Models <- factor(true_df$Models, levels=c("All models","Harley92,Ethier04&Xiao21",
                                              "vonC00,Dubois07&Busch18",
                                              "vonC00&Dubois07",
                                              "Tholen12","Yin04", "Busch18", "Xiao21"))

exp_vec <- c(expression(V[cmax]), expression(J[max]),
             expression(R[d]),expression(Gamma^"*"),
             expression(theta[J]),expression(alpha[J]),
             expression(K[C]),expression(K[O]),
             expression(g[m]),expression(V[TPU]),expression(alpha[G]),
             expression(Phi[IIm]),
             expression(f[Q]),expression(f[pseudo]),
             "h", expression(max[aG]),expression(max[aS]),
             expression(N[max]),expression(g[ch]),expression(g[wp]),
             "s",expression(Phi[IILL]))

levels(true_df$variable) <- exp_vec

g1 <- ggplot(true_df, aes(x = variable, y = value, fill = Models)) +
  geom_violin(trim = FALSE, alpha = 0.4, scale = "width") +
  labs(x = "", y = "Sampled values") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),  # rotate labels if many
    strip.text = element_text(size = 10)
  ) +
  facet_wrap(.~variable, scales="free",
             labeller = labeller(variable = label_parsed ),
             )

g1
#############################
load(file="results/Rdata/fig2_simulated_curves.RData")

# keep <- which(grepl("rep1_", names(simulated_curves), fixed = TRUE))
keep2 <- which(grepl("satQ", names(simulated_curves), fixed = TRUE))

# keep <- intersect(keep,keep2)
curves <- simulated_curves[keep2]


names(curves) <- sub("^rep[0-9]+_", "", names(curves))
names(curves) <- sub("_ACi\\(satQ\\)$", "", names(curves))

synth_curves <- map_dfr(curves, ~as.data.frame(.x), .id = "Model")

summary_df <- synth_curves %>%
  group_by(Model, Ci) %>%
  summarise(
    Median_A =round(median(An, na.rm = TRUE),2),
    Std_A = round(sd(An, na.rm = TRUE),3),
    lower_A = quantile(An, 0.05, na.rm = TRUE),
    upper_A = quantile(An, 0.95, na.rm = TRUE),
    .groups = 'drop'
  )

summary_df$Model[summary_df$Model=="Caemmerer00"] <- "vonC00"
# summary_df$Model[summary_df$Model=="Busch18"] <- "Busch17"
summary_df$Model[summary_df$Model=="Xiao20"] <- "Xiao21"

curves_name <- c("FvCB80", "Harley92","vonC00",
                 "Ethier04", "Yin04", "Dubois07",
                 "Tholen12", "Busch18", "Xiao21" )

summary_df$Model <- factor(summary_df$Model, levels = curves_name)

g2 <- ggplot(summary_df, aes(x = Ci)) +

  geom_errorbar(aes(ymin = lower_A, ymax = upper_A), width = 1, color = "orange") +
  geom_point(aes(y = Median_A),shape = 1, size = 2.5, color = "purple") +

  facet_wrap(~Model,ncol=3,dir = "h") +
  labs(
    x = expression(C[i]~"(µmol mol"^{-1}*")"),
    y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
    color = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10)) +
  coord_cartesian(ylim = c(-5,60))

ggarrange(g1,g2,labels = c('a', 'b'), ncol=1)


ggsave(filename = paste0("results/Figures/SFig1_sampled_params_and_simulation.png"),width = 7, height = 10)
