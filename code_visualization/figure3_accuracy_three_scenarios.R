library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(TeachingDemos)
setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")


# start <- T
# for (rep_x in 1:100){
#   file <- paste0("results/Rdata/fig2_samples_KCO/fig2_identifiability_KCO_i",rep_x,".RData")
#   if (file.exists(file)){
#     load(file)
#     for (i in 1:length(res_list)){
#       temp <- res_list[[i]]
#       ind <- which.max(temp$lp__)
#       map_vec <- temp[ind,]
#       rep <- paste0("rep",rep_x,"_")
#       rmind <- which(colnames(temp)%in%c("lp__",".draw"))
#       map_vec <- map_vec[,-rmind]
#
#       if (start){
#         MAP_df <- tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#           dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x)
#         start <- F
#       }else{
#         MAP_df <- bind_rows(MAP_df,tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#                               dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x))
#       }
#     }
#   }
# }
#
# save(MAP_df, file="results/Rdata/data_visualization/fig3_MAP_KCO.RData")



load(file="results/Rdata/data_visualization/fig3_MAP_KCO.RData")

load(file="results/Rdata/fig2_sampled_params.RData")


# 1) Extract all true parameter vectors (length 100 each)
true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")
# iter is 1..100, columns are params (Vcmax, Jmax, ...)

param_name <- c("V_cmax","J_max","R_d","gamma_star")
cases <-unique(MAP_df$scena_model)


cor_mat <- matrix(NA,length(cases),length(param_name))
colnames(cor_mat) <- param_name
rownames(cor_mat) <- cases


for (param in param_name){
  for (scena in cases){
    true_vec <- true_params_df[,param]
    est_vec <- MAP_df$value[MAP_df$scena_model==scena & MAP_df$param==param]
    if (length(est_vec)>0){
      cor_mat[scena, param] <- cor(true_vec,est_vec)
    }
  }
}


library(reshape2)

# convert correlation matrix to long format
df1 <- melt(cor_mat, varnames = c("model_scen", "param"), value.name = "correlation")
df1$estimation <- "Estimating all parameters"

##################################################################
# start <- T
# for (rep_x in 1:100){
#   file <- paste0("results/Rdata/fig3_samples_true_fixed/fig3_core_estimate_true_fixed_i",rep_x,".RData")
#   if (file.exists(file)){
#     load(file)
#     for (i in 1:length(res_list)){
#       temp <- res_list[[i]]
#       ind <- which.max(temp$lp__)
#       map_vec <- temp[ind,]
#       rep <- paste0("rep",rep_x,"_")
#       rmind <- which(colnames(temp)%in%c("lp__",".draw"))
#       map_vec <- map_vec[,-rmind]
#
#       if (start){
#         MAP_df <- tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#           dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x)
#         start <- F
#       }else{
#         MAP_df <- bind_rows(MAP_df,tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#                               dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x))
#       }
#     }
#   }
# }
#
# save(MAP_df, file="results/Rdata/data_visualization/fig3_MAP_true_fixed.RData")


load(file="results/Rdata/data_visualization/fig3_MAP_true_fixed.RData")

load(file="results/Rdata/fig2_sampled_params.RData")

true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")

param_name <- c("V_cmax","J_max","R_d","gamma_star")

cases <-unique(MAP_df$scena_model)


cor_mat <- matrix(NA,length(cases),length(param_name))
colnames(cor_mat) <- param_name
rownames(cor_mat) <- cases


for (param in param_name){
  for (scena in cases){
    true_vec <- true_params_df[,param]
    est_vec <- MAP_df$value[MAP_df$scena_model==scena & MAP_df$param==param]
    if (length(est_vec)>0){
      cor_mat[scena, param] <- cor(true_vec,est_vec)
    }
  }
}


df2 <- melt(cor_mat, varnames = c("model_scen", "param"), value.name = "correlation")
df2$estimation <- "Estimating core parameters\nnon-core fixed with true value"
##################################################################
# start <- T
# for (rep_x in 1:100){
#   file <- paste0("results/Rdata/fig3_samples_median_fixed/fig3_core_estimate_median_fixed_i",rep_x,".RData")
#   if (file.exists(file)){
#     load(file)
#     for (i in 1:length(res_list)){
#       temp <- res_list[[i]]
#       ind <- which.max(temp$lp__)
#       map_vec <- temp[ind,]
#       rep <- paste0("rep",rep_x,"_")
#       rmind <- which(colnames(temp)%in%c("lp__",".draw"))
#       map_vec <- map_vec[,-rmind]
#
#       if (start){
#         MAP_df <- tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#           dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x)
#         start <- F
#       }else{
#         MAP_df <- bind_rows(MAP_df,tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#                               dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x))
#       }
#     }
#   }
# }
#
# save(MAP_df, file="results/Rdata/data_visualization/fig3_MAP_median_fixed.RData")


load(file="results/Rdata/data_visualization/fig3_MAP_median_fixed.RData")

load(file="results/Rdata/fig2_sampled_params.RData")

true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")

param_name <- c("V_cmax","J_max","R_d","gamma_star")

cases <-unique(MAP_df$scena_model)


cor_mat <- matrix(NA,length(cases),length(param_name))
colnames(cor_mat) <- param_name
rownames(cor_mat) <- cases


for (param in param_name){
  for (scena in cases){
    true_vec <- true_params_df[,param]
    est_vec <- MAP_df$value[MAP_df$scena_model==scena & MAP_df$param==param]
    if (length(est_vec)>0){
      cor_mat[scena, param] <- cor(true_vec,est_vec)
    }
  }
}


df3 <- melt(cor_mat, varnames = c("model_scen", "param"), value.name = "correlation")
df3$estimation <- "Estimating core parameters\nnon-core fixed with median value"


df <- rbind(df1,df2)
df <- rbind(df,df3)

df <- df[df$model_scen!="Xiao20_ACi(different Q) & CF",]
df <- df %>%
  separate(model_scen, into = c("model", "datatype"), sep = "_", extra = "merge")

df$model[df$model=="Caemmerer00"] <- "vonC.00"
df$model[df$model=="Busch18"] <- "Busch17"
df$model[df$model=="Xiao20"] <- "Xiao21"
df$datatype[df$datatype=="ACi(satQ)"] <- "Single A-Ci"
df$datatype[df$datatype=="ACi & AQ"] <- "A-Ci & A-Q"
df$datatype[df$datatype=="ACi(different Q)"] <- "Multiple A-Ci"


df$param <- factor(df$param,levels=param_name)

newmodels <-  c("FvCB80", "Harley92","vonC.00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch17","Xiao21")

df$model <- factor(df$model,levels=newmodels)

datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")
df$datatype <- factor(df$datatype, levels=datatypes)


check <- df %>% group_by(param,datatype) %>%
  summarise(mean_correlation = mean(correlation, na.rm = TRUE))


exp_vec <- c(expression(V[cmax]),
             expression(J[max]),
             expression(R[d]),
             expression(Gamma^"*"))

levels(df$param)<- exp_vec


# ggplot(df, aes(x = model, y = correlation, fill = estimation)) +
#   geom_col(position = position_dodge(width = 0.85), width = 0.75) +
#   facet_grid(param ~ datatype, scales = "free_x", space = "free_x") +
#   coord_cartesian(ylim = c(0, 1)) +
#   labs(x = NULL, y = "Correlation (estimated vs true)", fill = "Estimation") +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.major.x = element_blank(),
#     legend.position = "bottom"
#   )


ggplot(df, aes(x = model, y = correlation, color = estimation, shape=estimation)) +
  geom_point(size = 2) +
  geom_line(aes(group = model), linewidth = 0.4, alpha = 0.5) +
  facet_grid(param ~ datatype,
             switch = "y",
             labeller = labeller(param = label_parsed)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Correlation", color = "", shape= "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size=9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(filename = paste0("results/Figures/Fig3_accuracy_new.png"),width = 7, height = 6)

# g2
# library("ggpubr")
# ggarrange(g1,g2,labels = c('a', 'b'), ncol=1) #, heights = c(0.8,1)
#
# ggsave(filename = paste0("results/Figures/Fig3_accuracy.png"),width = 7, height = 6)
#
# ggarrange(g1_t,g2_t,labels = c('a', 'b'), ncol=1) #, heights = c(0.8,1)
#
# ggsave(filename = paste0("results/Figures/Fig3_accuracy_text.png"),width = 7, height = 6)
#
