library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(ggpubr)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")

ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
ref <- read.csv("data/2022_barley_reference.csv")

AQ <- AQ[AQ$Flag_removal!="x",]
ACi$genotype <- ref$Accession[match(ACi$Plot,ref$PlotID)]
AQ$genotype <- ref$Accession[match(AQ$Plot,ref$PlotID)]

genotype_id <- intersect(ACi$genotype, AQ$genotype)


datatypes <- c("ACi", "ACi & AQ")
com_pars <- c("V_cmax", "J_max", "R_d", "gamma_star","theta_J","alpha_J","K_CO")

models <- c("FvCB80", "Harley92", "vonC.00", "Ethier04", "Yin04",
            "Dubois07", "Tholen12", "Busch17", "Xiao21" )
models0 <- c("FvCB80", "Harley92", "Caemmerer00", "Ethier04", "Yin04",
            "Dubois07", "Tholen12", "Busch18", "Xiao20" )

mse      <- function(y, yhat) mean((y - yhat)^2, na.rm = TRUE)
rmse     <- function(y, yhat) sqrt(mse(y, yhat))

my_cv <- function(x){
  sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
}

df_dist <- list()
df_var <- list()
cor_list <- list()

for (t in 1:2){
  dist_list_norm <- list()
  dist_list <- list()
  for (rep_x in 1:length(genotype_id)){

    param_matrix <- matrix(NA, nrow=length(com_pars), ncol=9)
    rownames(param_matrix) <- com_pars
    colnames(param_matrix) <- models

    file <- paste0("results/Rdata/fig4_genotype_KCO_B/fig4_genotype_fit_i",rep_x,".RData")
    load(file)

    if (t==1){
      idx <- grep("_ACi$", names(res_list))
    }else{
      idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
    }

    res_list <- res_list[idx]
    res_models <- sub("^[^_]+_([^_]+)_.*$", "\\1", names(res_list))

    for (i in 1:length(res_list)){
      temp <- res_list[[i]]

      if (nrow(temp)!=1 & ncol(temp)!=10){
        ind <- which.max(temp$lp__)
        map_vec <- temp[ind,]

        model <- models[match(res_models[i], models0)]

        keep <- intersect(com_pars, colnames(map_vec))
        param_matrix[match(keep, com_pars),model] <- as.numeric(map_vec[1, keep])
      }

    }

    norm_matrix <- param_matrix #scale(param_matrix)
    dist_list[[rep_x]] <- dist(t(norm_matrix), method = "euclidean")
    dist_list_norm[[rep_x]] <- dist(t(scale(norm_matrix)), method = "euclidean")
  }

  arr_norm <- simplify2array(lapply(dist_list_norm, as.matrix))
  mean_mat <- apply(arr_norm, c(1, 2), mean, na.rm = TRUE)
  cor_list[[t]] <- mean_mat
  mean_mat[upper.tri(mean_mat, diag = TRUE)] <- NA

  arr <- simplify2array(lapply(dist_list, as.matrix))
  var_mat <- apply(arr, c(1, 2), my_cv)
  var_mat[upper.tri(var_mat, diag = TRUE)] <- NA

  rownames(mean_mat) <- models
  colnames(mean_mat) <- models

  rownames(var_mat) <- models
  colnames(var_mat) <- models

  df_dist[[t]] <- melt(mean_mat)
  colnames(df_dist[[t]]) <- c("Model1", "Model2","Distance")

  df_dist[[t]]$CV <- melt(var_mat)$value
  df_dist[[t]]$Data <- datatypes[t]
}

df_dist0 <- df_dist
df_dist <- rbind(df_dist[[1]],df_dist[[2]])


library(pheatmap)

g1 <- pheatmap(cor_list[[1]],
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         color = colorRampPalette(c("#77DD77", "#FFF176", "#FFB347"))(50),
         display_numbers = TRUE,
         number_format = "%.2f",
         angle_col =  "45")

g2 <- pheatmap(cor_list[[2]],
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "average",
               color = colorRampPalette(c("#77DD77", "#FFF176", "#FFB347"))(50),
               display_numbers = TRUE,
               number_format = "%.2f",
               angle_col =  "45")

library(ggplotify)
library(ggpubr)

p1 <- as.ggplot(g1$gtable)
p2 <- as.ggplot(g2$gtable)

ggarrange(p1, p2, labels = c("a", "b"), ncol = 1)

ggsave(filename = paste0("results/Figures/SFig9_model_dendrogram.png"),width = 6, height = 7)

