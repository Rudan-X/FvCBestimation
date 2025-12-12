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

common_par <-c("V_cmax","J_max","R_d","gamma_star") #,"theta_J","alpha_J", "K_C", "K_O"


newmodels <-  c("FvCB80", "Harley92","vonCaemmerer00",
              "Ethier04", "Yin04","Dubois07",
              "Tholen12", "Busch17","Xiao21")


models <-  c("FvCB80", "Harley92","Caemmerer00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch18","Xiao20")

# start <- T
# for (rep_x in 1:100){
#   file <- paste0("results/Rdata/fig2_samples/fig1_identifiability_i",rep_x,".RData")
#   load(file)
#   true_par <- as.data.frame(sampled_params[[rep_x]])
#
#   for (i in 1:length(res_list)){
#     temp_df <- res_list[[i]]
#
#     rmind <- which(colnames(temp_df)%in%c("lp__",".draw"))
#     temp_df <- temp_df[,-rmind]
#
#     model_scena <- names(res_list)[i]
#     check <- strsplit(model_scena,"_")
#     model_x <- check[[1]][2]
#     # model_x <- newmodels[match(model_x,models)]
#     dat_x <- check[[1]][3]
#
#
#     for (t in 1: ncol(temp_df)){
#       CIs <- hpd_interval(as.numeric(unlist(temp_df[,t])),prob=0.5)
#       var_x <- colnames(temp_df)[t]
#       range_x <- bounds[[var_x]][2]-bounds[[var_x]][1]
#       CI_range <- CIs[2] - CIs[1]
#
#       temp <- data.frame(variable=var_x,model=model_x,data=dat_x, range_perc=CI_range/range_x*100)
#       if (start){
#         df <- temp
#         start <- F
#       }else{
#         df <- rbind(df,temp)
#       }
#     }
#   }
# }
#
# save(df, file = "results/Rdata/data_visualization/fig2_CIranges.RData")

load( file = "results/Rdata/data_visualization/fig2_CIranges.RData")
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

exp_vec <- c(expression(V[cmax]),
             expression(J[max]),
             expression(R[d]),
             expression(Gamma^"*"))
# ,
#              expression(theta[J]),
#              expression(alpha[J]),
#              expression(K[C]),
#              expression(K[O]))
levels(df$variable)<- exp_vec

cols <- c("#1B9E77", "#D95F02", "#7570B3")
ggplot(df,
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
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size=11),
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



ggsave(filename = paste0("results/Figures/Fig2.png"),width = 7, height = 6)
