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


param_name <- colnames(true_params_df)
param_name <- param_name[-1]

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
df <- melt(cor_mat, varnames = c("model_scen", "param"), value.name = "correlation")

df <- df[df$model_scen!="Xiao20_ACi(different Q) & CF",]
df <- df %>%
  separate(model_scen, into = c("model", "datatype"), sep = "_", extra = "merge")

df$model[df$model=="Caemmerer00"] <- "vonC.00"
df$model[df$model=="Busch18"] <- "Busch17"
df$model[df$model=="Xiao20"] <- "Xiao21"
df$datatype[df$datatype=="ACi(satQ)"] <- "Single A-Ci"
df$datatype[df$datatype=="ACi & AQ"] <- "A-Ci & A-Q"
df$datatype[df$datatype=="ACi(different Q)"] <- "Multiple A-Ci"


common_par <-c("V_cmax","J_max","R_d","gamma_star") # ,"theta_J","alpha_J", "K_CO"

common_cor <- df[df$param%in%common_par,]
common_cor$param <- factor(common_cor$param,levels=rev(common_par))

newmodels <-  c("FvCB80", "Harley92","vonC.00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch17","Xiao21")

common_cor$model <- factor(common_cor$model,levels=newmodels)

datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")
common_cor$datatype <- factor(common_cor$datatype,
                              levels=datatypes) #"GE & CF"


# common_cor$correlation[common_cor$param=="theta_J" & common_cor$model=="Harley92"] <- NA

check <- common_cor %>% group_by(param,datatype) %>%
  summarise(mean_correlation = mean(correlation, na.rm = TRUE))

mincor=0.19
maxcor=1
midp <- mincor+(maxcor-mincor)/2

exp_vec <- c(expression(V[cmax]),
             expression(J[max]),
             expression(R[d]),
             expression(Gamma^"*"))

# ,
# expression(theta[J]),
# expression(alpha[J]),
# expression(K[C]),
# expression(K[O])

g1 <- ggplot(common_cor, aes(x = model, y = param, fill = correlation, label = round(correlation, 2))) +
  geom_tile(color = "white") +
  # geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(mincor, maxcor), name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom") +
  facet_grid(.~datatype) +
  scale_y_discrete(labels = rev(exp_vec))

g1_t <- ggplot(common_cor, aes(x = model, y = param, fill = correlation, label = round(correlation, 2))) +
  geom_tile(color = "white") +
  geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(mincor, maxcor), name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom") +
  facet_grid(.~datatype) +
  scale_y_discrete(labels = rev(exp_vec))
g1

##################################################################
#
# MAP_list <- list()
# for (r in 1:5){
#   start <- T
#
#   finished <- c()
#   for (rep_x in 1:100){
#     print(paste0("R", r, ", rep_x: ", rep_x))
#     file <- paste0("results/Rdata/fig3_samples/fig3_core_estimate_i",rep_x,".RData")
#     if (file.exists(file)){
#       load(file)
#       finished <- c(finished,rep_x)
#       for (i in 1:length(res_list)){
#         temp <- res_list[[i]][[r]]
#         ind <- which.max(temp$lp__)
#         map_vec <- temp[ind,]
#         rep <- paste0("rep",rep_x,"_")
#         rmind <- which(colnames(temp)%in%c("lp__",".draw"))
#         map_vec <- map_vec[,-rmind]
#
#         if (start){
#           temp_df <- tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#             dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x)
#           start <- F
#         }else{
#           temp_df <- bind_rows(temp_df,tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#                                  dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x))
#         }
#       }
#     }
#   }
#   MAP_list[[r]] <- temp_df
# }
# save(MAP_list, file="results/Rdata/data_visualization/Fig3_MAPestimates.RData")

# MAP_list2 <- list()
# for (r in 1:5){
#   start <- T
#   finished <- c()
#   for (rep_x in 1:100){
#     print(paste0("R", r, ", rep_x: ", rep_x))
#     file <- paste0("results/Rdata/fig3_samples_v2/fig3_core_estimate_i",rep_x,".RData")
#     if (file.exists(file)){
#       load(file)
#       finished <- c(finished,rep_x)
#       for (i in 1:length(res_list)){
#         temp <- res_list[[i]][[r]]
#         ind <- which.max(temp$lp__)
#         map_vec <- temp[ind,]
#         rep <- paste0("rep",rep_x,"_")
#         rmind <- which(colnames(temp)%in%c("lp__",".draw"))
#         map_vec <- map_vec[,-rmind]
#
#         if (start){
#           temp_df <- tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#             dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x)
#           start <- F
#         }else{
#           temp_df <- bind_rows(temp_df,tidyr::pivot_longer(map_vec, everything(), names_to = "param", values_to = "value") |>
#                                  dplyr::mutate(scena_model = gsub(rep,"",names(res_list)[i]), rep = rep_x))
#         }
#       }
#     }
#   }
#   MAP_list2[[r]] <- temp_df
# }
# save(MAP_list2, file="results/Rdata/data_visualization/Fig3_MAPestimates6to10.RData")

load( file="results/Rdata/data_visualization/Fig3_MAPestimates.RData")
load( file="results/Rdata/data_visualization/Fig3_MAPestimates6to10.RData")


MAP_list <- c(MAP_list, MAP_list2)
load(file="results/Rdata/fig2_sampled_params.RData")
true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")
param_name <- c("V_cmax","J_max","R_d","gamma_star")

cor_list <- list()
for (r in 1:10){
  cases <-unique(MAP_list[[r]]$scena_model)

  cor_mat <- matrix(NA,length(cases),length(param_name))
  colnames(cor_mat) <- param_name
  rownames(cor_mat) <- cases

  for (param in param_name){
    for (scena in cases){
      true_vec <- true_params_df[,param]
      est_vec <- MAP_list[[r]]$value[MAP_list[[r]]$scena_model==scena & MAP_list[[r]]$param==param]
      if (length(est_vec)>0){
        cor_mat[scena, param] <- cor(true_vec,est_vec)
      }
    }
  }
  cor_list[[r]] <- cor_mat
}

library(abind)

cor_array <- abind::abind(cor_list, along = 3)

cor_avg <- apply(cor_array, c(1, 2), mean)


library(reshape2)


df <- melt(cor_avg, varnames = c("model_scen", "param"), value.name = "correlation")
df <- df[df$model_scen!="Xiao20_ACi(different Q) & CF",]
df <- df %>%
  separate(model_scen, into = c("model", "datatype"), sep = "_", extra = "merge")

df$model[df$model=="Caemmerer00"] <- "vonC.00"
df$model[df$model=="Busch18"] <- "Busch17"
df$model[df$model=="Xiao20"] <- "Xiao21"
df$datatype[df$datatype=="ACi(satQ)"] <- "Single A-Ci"
df$datatype[df$datatype=="ACi & AQ"] <- "A-Ci & A-Q"
df$datatype[df$datatype=="ACi(different Q)"] <- "Multiple A-Ci"


exp_vec <- c(expression(V[cmax]),
             expression(J[max]),
             expression(R[d]),
             expression(Gamma^"*"))


df$param <- factor(df$param,levels=rev(param_name))
newmodels <-  c("FvCB80", "Harley92","vonC.00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch17","Xiao21")
df$model <- factor(df$model,levels=newmodels)
datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")
df$datatype <- factor(df$datatype,levels=datatypes)

mincor=0.19
maxcor=1
midp <- mincor+(maxcor-mincor)/2


g2 <- ggplot(df, aes(x = model, y = param, fill = correlation, label = round(correlation, 2))) +
  geom_tile(color = "white") +
  # geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(mincor, maxcor), name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "none") +
  facet_grid(.~datatype,  drop = TRUE) +
  scale_y_discrete(labels = rev(exp_vec))

g2_t <- ggplot(df, aes(x = model, y = param, fill = correlation, label = round(correlation, 2))) +
  geom_tile(color = "white") +
  geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(mincor, maxcor), name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "none") +
  facet_grid(.~datatype,  drop = TRUE) +
  scale_y_discrete(labels = rev(exp_vec))


df_summary <- df %>%
  group_by(param, datatype) %>%
  summarise(mean_correlation = mean(correlation, na.rm = TRUE)) %>%
  ungroup()

df_summary

g2
library("ggpubr")
ggarrange(g1,g2,labels = c('a', 'b'), ncol=1) #, heights = c(0.8,1)

ggsave(filename = paste0("results/Figures/Fig3_accuracy.png"),width = 7, height = 6)

ggarrange(g1_t,g2_t,labels = c('a', 'b'), ncol=1) #, heights = c(0.8,1)

ggsave(filename = paste0("results/Figures/Fig3_accuracy_text.png"),width = 7, height = 6)

