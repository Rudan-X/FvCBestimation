library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(TeachingDemos)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")
load(file="results/Rdata/fig2_sampled_params.RData")
load(file="results/Rdata/data_visualization/Sfig_TPU_ind.RData")


# hpd_interval <- function(x, prob = 0.95) {
#   x <- sort(as.numeric(x))
#   n <- length(x)
#   k <- floor(prob * n)
#
#   if (k < 1L || n <= k) stop("Not enough samples")
#
#   # sliding window of size k+1 (covers 'prob' mass)
#   widths <- x[(k+1):n] - x[1:(n-k)]
#   i_min  <- which.min(widths)
#
#   lower <- x[i_min]
#   upper <- x[i_min + k]
#
#   c(lower = lower, upper = upper)
# }
#
# bounds <- list(
#   V_cmax     = c(50, 200),
#   J_max      = c(100, 250),
#   R_d        = c(0.5, 5),
#   gamma_star = c(10, 35),
#   V_tpu      = c(8, 20),
#   gm         = c(0.1, 0.7),
#   alpha_J    = c(0.3, 0.5),
#   theta_J    = c(0.5, 0.9),
#   K_C        = c(200, 300),
#   K_O        = c(150, 300),
#   alpha_G    = c(0.3, 0.8),
#   phi2m      = c(0.5, 0.9),
#   f_Q        = c(0.5, 1),
#   f_pseudo   = c(0, 0.3),
#   h          = c(3, 14/3),
#   max_aG     = c(0.01, 0.2),
#   max_aS     = c(0.01, 0.5),
#   N_max      = c(0.1, 2),
#   g_ch       = c(0.05, 0.7),
#   g_wp       = c(0.1, 2),
#   s          = c(0.3, 0.6),
#   Phi2LL     = c(0.5, 1)
# )
#
#
# start <- T
# for (rep_x in tpu_ind){
#   file <- paste0("../FvCBestimation_not_uploaded_true/results/Rdata/fig2_samples/fig1_identifiability_i",rep_x,".RData")
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
#         df_tpu <- temp
#         start <- F
#       }else{
#         df_tpu <- rbind(df_tpu,temp)
#       }
#     }
#   }
# }
#
# save(df_tpu, file = "results/Rdata/data_visualization/fig2_CIranges_tpu_only.RData")
#
#
# start <- T
# for (rep_x in setdiff(1:100,tpu_ind)){
#   file <- paste0("../FvCBestimation_not_uploaded_true/results/Rdata/fig2_samples/fig1_identifiability_i",rep_x,".RData")
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
#         df_notpu <- temp
#         start <- F
#       }else{
#         df_notpu <- rbind(df_notpu,temp)
#       }
#     }
#   }
# }
#
# save(df_notpu, file = "results/Rdata/data_visualization/fig2_CIranges_notpu.RData")
###################################################################################
newmodels <-  c("FvCB80", "Harley92","vonC00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch18","Xiao21")


models <-  c("FvCB80", "Harley92","Caemmerer00",
             "Ethier04", "Yin04","Dubois07",
             "Tholen12", "Busch18","Xiao20")


#########################################
load(file = "results/Rdata/data_visualization/fig2_CIranges_notpu.RData")
df <- df_notpu

df <- df[df$data=="ACi(satQ)",]
df$data[df$data=="ACi(satQ)"] <- "Single A-Ci"



ind <- which(df$variable=="theta_J" & df$model=="Harley92")
df <- df[-ind,]


df$model <- newmodels[match(df$model,models)]
df$model <- factor(df$model, levels = newmodels)



specific_par <- c("V_tpu","alpha_G","max_aG","max_aS","N_max")

df <- df[df$variable%in%specific_par,]
df$variable <- factor(df$variable, levels = specific_par)

df$model <- factor(df$model, levels = newmodels)


var_map <- c(
  V_tpu    = "V[TPU]",
  alpha_G  = "alpha[G]",
  max_aG   = "max[aG]",
  max_aS   = "max[aS]",
  N_max    = "N[max]"
)


combos <- tibble::tribble(
  ~model,     ~variable,
  "vonC00","V_tpu",
  "vonC00","alpha_G",

  "Dubois07","V_tpu",
  "Dubois07","alpha_G" ,

  "Busch18","V_tpu",
  "Busch18","max_aG",
  "Busch18","max_aS",
  "Busch18","N_max"
)

df1 <- semi_join(df, combos, by = c("model","variable")) %>%
  mutate(
    # make a model • variable label where variable is plotmath
    variable_lab = dplyr::recode(variable, !!!var_map),
    facet_lab = paste(model,"~", variable_lab),  # plain text + plotmath
    facet_lab = factor(
      facet_lab,
      levels = paste(combos$model,"~", dplyr::recode(combos$variable, !!!var_map))
    )
  )

df1$data <- factor(df1$data, levels = datatypes)
df1$model <- factor(df1$model, levels = newmodels)
df1$type <- "no TPU limitation present"

#########################################
load(file = "results/Rdata/data_visualization/fig2_CIranges_tpu_only.RData")
df <- df_tpu

df <- df[df$data=="ACi(satQ)",]
df$data[df$data=="ACi(satQ)"] <- "Single A-Ci"

ind <- which(df$variable=="theta_J" & df$model=="Harley92")
df <- df[-ind,]


df$model <- newmodels[match(df$model,models)]
df$model <- factor(df$model, levels = newmodels)




specific_par <- c("V_tpu","alpha_G","max_aG","max_aS","N_max")

df <- df[df$variable%in%specific_par,]
df$variable <- factor(df$variable, levels = specific_par)

df$model <- factor(df$model, levels = newmodels)


var_map <- c(
  V_tpu    = "V[TPU]",
  alpha_G  = "alpha[G]",
  max_aG   = "max[aG]",
  max_aS   = "max[aS]",
  N_max    = "N[max]"
)


combos <- tibble::tribble(
  ~model,     ~variable,
  "vonC00","V_tpu",
  "vonC00","alpha_G",

  "Dubois07","V_tpu",
  "Dubois07","alpha_G" ,

  "Busch18","V_tpu",
  "Busch18","max_aG",
  "Busch18","max_aS",
  "Busch18","N_max"
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
df2$type <- "TPU limitation present"


df <- rbind(df1,df2)

cols <- c("#F4A6A1", "#9B7EDC")

g2 <- ggplot(df, aes(y= range_perc, x = data, fill = type, colour = type)) +
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


library("ggpubr")
ggarrange(g1,g2,labels = c('a', 'b'), ncol=1, heights = c(0.8,1))

ggsave(filename = paste0("results/Figures/SFig4_TPUeffect.png.png"),width = 7, height = 8)



###############


ggplot(df, aes(y= range_perc, x = data)) +
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
  labs(x = "", y = "Relative CI width (%)") +
  coord_cartesian(ylim = c(0, 65))+
  scale_y_continuous(breaks = seq(0, 60, 20))  # specify breaks

