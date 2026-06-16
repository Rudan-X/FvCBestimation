library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation_not_uploaded_true/")

ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
ACi$ID <- paste0(ACi$Plot,"_",ACi$Repeat)
AQ$ID <- paste0(AQ$Plot,"_",ACi$Repeat)

plant_id <- intersect(ACi$ID, AQ$ID)

common_plot <- intersect(ACi$Plot, AQ$Plot)
ACi <- ACi[ACi$Plot%in%common_plot,]
AQ <- AQ[AQ$Plot%in%common_plot,]

datatypes <- c("ACi", "ACi & AQ")
com_pars <- c("V_cmax", "J_max", "R_d", "gamma_star","theta_J","alpha_J","K_CO",
              "gm","V_tpu","alpha_G", "phi2m", "f_Q","f_pseudo", "h","max_aG","max_aS","N_max",
              "g_ch",  "g_wp","s","Phi2LL" )


models <- c("FvCB80", "Harley92", "vonC00", "Ethier04", "Yin04",
            "Dubois07", "Tholen12", "Busch18", "Xiao21" )

models0 <- c("FvCB80", "Harley92", "Caemmerer00", "Ethier04", "Yin04",
            "Dubois07", "Tholen12", "Busch18", "Xiao20" )

plot_list <- list()
plant_list <- list()
plant_var_list <- list()


for (rep_x in 1:length(common_plot)){ #
  file <- paste0("results/Rdata/fig4_plot_KCO_B/fig4_plot_KCO_i",rep_x,".RData")
  load(file)

  idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
  res_list <- res_list[idx]

  res_models <- sub("^[^_]+_([^_]+)_.*$", "\\1", names(res_list))

  param_matrix <- matrix(NA, nrow=length(com_pars), ncol=9)
  rownames(param_matrix) <- com_pars

  for (i in 1:length(res_list)){
    temp <- res_list[[i]]

    if (nrow(temp)==1 & ncol(temp)==10){
      param_matrix[, match(res_models[i], models0)] <- NA
    }else{
      ind <- which.max(temp$lp__)
      map_vec <- temp[ind,]

      keep <- intersect(com_pars, colnames(map_vec))        # common params
      param_matrix[match(keep, com_pars), match(res_models[i], models0)] <- as.numeric(map_vec[1, keep])
    }

    # if ((rep_x %in%c(35,  74,  96, 100, 116)) & (i==1)){
    #   param_matrix[match(keep, com_pars), 1] <- as.numeric(map_vec[1, keep])
    # }else{
    #
    # }

  }
  plot_list[[rep_x]] <- param_matrix
  target_plot <- as.character(common_plot[rep_x])  # e.g. "1001"
  pind <- grep(paste0("^", target_plot, "_"), plant_id)

  if (length(pind)!=0){
    arr <- array(NA, c(length(com_pars),9,length(pind)))
    for (id in 1:length(pind)){
      file <- paste0("results/Rdata/fig4_plant_KCO_B/fig4_plant_fit_i",pind[id],".RData")
      load(file)
      idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
      res_list <- res_list[idx]

      res_models <- sub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", names(res_list))
      for (i in 1:length(res_list)){
        temp <- res_list[[i]]
        if (nrow(temp)==1 & ncol(temp)==10){
          arr[, match(res_models[i], models0), id] <- NA
        }else{
          ind <- which.max(temp$lp__)
          map_vec <- temp[ind,]

          keep <- intersect(com_pars, colnames(map_vec))        # common params
          arr[match(keep, com_pars), match(res_models[i], models0), id] <- as.numeric(map_vec[1, keep])
        }

      }
    }
    plant_list[[rep_x]] <- apply(arr, c(1,2), mean, na.rm= TRUE)
    plant_var_list[[rep_x]] <- apply(arr, c(1, 2), var, na.rm = TRUE)
  }else{
    plant_list[[rep_x]] <- matrix(NA, nrow=length(com_pars),ncol=9)
    plant_var_list[[rep_x]] <- matrix(NA, nrow=length(com_pars),ncol=9)
  }
}

plant_matrix <- simplify2array(lapply(plant_list, as.matrix))
plot_matrix <- simplify2array(lapply(plot_list, as.matrix))

cor_mat <- matrix(NA, length(com_pars), length(models))

for (p in 1:length(com_pars)){
  for (m in 1:length(models)){
    cor_mat[p,m] <- cor(plant_matrix[p,m,],plot_matrix[p,m,], use = "pairwise.complete.obs")
  }
}

rownames(cor_mat) <- com_pars
colnames(cor_mat) <- models

corr_df1 <- melt(cor_mat)
colnames(corr_df1) <- c("parameter","model","correlation")
corr_df1$compare <- "Plot- vs plant-level"


##########################################

ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
ref <- read.csv("data/2022_barley_reference.csv")

AQ <- AQ[AQ$Flag_removal!="x",]
ACi$genotype <- ref$Accession[match(ACi$Plot,ref$PlotID)]
AQ$genotype <- ref$Accession[match(AQ$Plot,ref$PlotID)]
ACi$ID <- paste0(ACi$Plot,"_",ACi$Repeat)
AQ$ID <- paste0(AQ$Plot,"_",ACi$Repeat)


genotype_id <- intersect(ACi$genotype, AQ$genotype)

ACi <- ACi[ACi$genotype%in%genotype_id,]
AQ <- AQ[AQ$genotype%in%genotype_id,]


datatypes <- c("ACi", "ACi & AQ")
com_pars <- c("V_cmax", "J_max", "R_d", "gamma_star","theta_J","alpha_J","K_CO",
              "gm","V_tpu","alpha_G", "phi2m", "f_Q","f_pseudo", "h",
              "g_ch",  "g_wp","max_aG","max_aS","N_max","s","Phi2LL" )


geno_list <- list()
plant_list <- list()
plant_var_list <- list()


for (rep_x in 1:length(genotype_id)){ #
  file <- paste0("results/Rdata/fig4_genotype_KCO_B/fig4_genotype_fit_i",rep_x,".RData")
  load(file)
  idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
  res_list <- res_list[idx]
  res_models <- sub("^[^_]+_([^_]+)_.*$", "\\1", names(res_list))

  param_matrix <- matrix(NA, nrow=length(com_pars), ncol=9)
  rownames(param_matrix) <- com_pars

  for (i in 1:length(res_list)){
    temp <- res_list[[i]]
    if (nrow(temp)==1 & ncol(temp)==10){
      param_matrix[, match(res_models[i], models0)] <- NA
    }else{
      ind <- which.max(temp$lp__)
      map_vec <- temp[ind,]

      keep <- intersect(com_pars, colnames(map_vec))        # common params
      param_matrix[match(keep, com_pars), match(res_models[i], models0)] <- as.numeric(map_vec[1, keep])
    }
  }
  geno_list[[rep_x]] <- param_matrix
  plots <- ref$PlotID[ref$Accession==genotype_id[rep_x]]
  target_plot <- as.character(plots)  # e.g. "1001"

  pind <- c()
  for (k in target_plot){
    pind <- c(pind,grep(paste0("^", k, "_"), plant_id))
  }

  if (length(pind)!=0){
    arr <- array(NA, c(length(com_pars),9,length(pind)))
    for (id in 1:length(pind)){
      file <- paste0("results/Rdata/fig4_plant_KCO_B/fig4_plant_fit_i",pind[id],".RData")
      load(file)
      idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
      res_list <- res_list[idx]
      res_models <- sub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", names(res_list))
      for (i in 1:length(res_list)){
        temp <- res_list[[i]]
        if (nrow(temp)==1 & ncol(temp)==10){
          arr[, match(res_models[i], models0), id] <- NA
        }else{
          ind <- which.max(temp$lp__)
          map_vec <- temp[ind,]

          keep <- intersect(com_pars, colnames(map_vec))        # common params
          arr[match(keep, com_pars), match(res_models[i], models0), id] <- as.numeric(map_vec[1, keep])
        }
      }
    }
    plant_list[[rep_x]] <- apply(arr, c(1,2), mean, na.rm= TRUE)
    plant_var_list[[rep_x]] <- apply(arr, c(1, 2), var, na.rm = TRUE)
  }else{
    plant_list[[rep_x]] <- matrix(NA, nrow=length(com_pars),ncol=9)
    plant_var_list[[rep_x]] <- matrix(NA, nrow=length(com_pars),ncol=9)
  }
}

plant_matrix <- simplify2array(lapply(plant_list, as.matrix))
geno_matrix <- simplify2array(lapply(geno_list, as.matrix))

cor_mat <- matrix(NA, length(com_pars), length(models))
for (p in 1:length(com_pars)){
  for (m in 1:length(models)){
    cor_mat[p,m] <- cor(plant_matrix[p,m,],geno_matrix[p,m,], use = "pairwise.complete.obs")
  }
}

rownames(cor_mat) <- com_pars
colnames(cor_mat) <- models

corr_df2 <- melt(cor_mat)
colnames(corr_df2) <- c("parameter","model","correlation")
corr_df2$compare <- "Genotype- vs plant-level"

exp_vec <- c(expression(italic(V[cmax])), expression(italic(J[max])),
             expression(italic(R[d])),expression(italic("Γ*")),
             expression(italic(theta[J])),
             expression(italic(alpha[J])),
             expression(italic(K[CO])),
             expression(italic(g[m])),expression(italic(V[TPU])),expression(italic(alpha[G])),
             expression(italic(Phi[IIm])),expression(italic(f[Q])),expression(italic(f[pseudo])),expression(italic("h")),
             expression(italic(g[ch])),expression(italic(g[wp])),
             expression(italic(max[aG])),expression(italic(max[aS])),expression(italic(N[max])),
             expression(italic("s")),expression(italic(Phi[IILL])))


corr_df <- rbind(corr_df1,corr_df2)

corr_df$parameter <- factor(corr_df$parameter,levels=rev(com_pars))
# corr_df$compare <- factor(corr_df$compre,levels=c("Plot- vs plant-level","Genotype- vs plant-level"))

#####################
mincor=max(corr_df$correlation, na.rm = T)
maxcor=min(corr_df$correlation, na.rm = T)
midp <- mincor+(maxcor-mincor)/2

ggplot(corr_df, aes(x = model, y = parameter, fill = correlation, label = round(correlation, 2))) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(mincor, maxcor), name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(face="bold")) +
  facet_grid(.~compare) +
  scale_y_discrete(labels = rev(exp_vec))

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")
ggsave(filename = paste0("results/Figures/SFig9_plant_vs_plot.png"),width = 6, height = 7)



ggplot(corr_df, aes(x = model, y = parameter, fill = correlation, label = round(correlation, 2))) +
  geom_tile(color = "white") +
  # geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(mincor, maxcor), name = "Correlation") +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(face="bold")) +
  facet_grid(.~compare) +
  scale_y_discrete(labels = rev(exp_vec))


ggsave(filename = paste0("results/Figures/SFig8_plant_vs_plot_text.png"),width = 6, height = 7)





