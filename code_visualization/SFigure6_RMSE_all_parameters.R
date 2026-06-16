library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")
mse      <- function(y, yhat) mean((y - yhat)^2, na.rm = TRUE)
rmse     <- function(y, yhat) sqrt(mse(y, yhat))


load(file="results/Rdata/data_visualization/fig3_MAP_KCO.RData")

load(file="results/Rdata/fig2_sampled_params.RData")
for (i in 1:100){
  sampled_params[[i]]$K_CO <- sampled_params[[i]]$K_C * (1 + 210 / sampled_params[[i]]$K_O)
}


# 1) Extract all true parameter vectors (length 100 each)
true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")
# iter is 1..100, columns are params (Vcmax, Jmax, ...)


param_name <- c("V_cmax","J_max","R_d","gamma_star","theta_J","alpha_J", "K_CO",
                "gm","V_tpu","alpha_G", "phi2m", "f_Q","f_pseudo", "h",
                "g_ch",  "g_wp","max_aG","max_aS","N_max","s","Phi2LL")

cases <-unique(MAP_df$scena_model)


RMSE_mat <- matrix(NA,length(cases),length(param_name))
colnames(RMSE_mat) <- param_name
rownames(RMSE_mat) <- cases


for (param in param_name){
  for (scena in cases){
    true_vec <- true_params_df[,param]
    est_vec <- MAP_df$value[MAP_df$scena_model==scena & MAP_df$param==param]
    if (length(est_vec)>0){
      RMSE_mat[scena, param] <- sqrt(mean(((est_vec - true_vec) / true_vec)^2))# rmse(true_vec, est_vec)
    }
  }
}

library(reshape2)

# convert correlation matrix to long format
df <- melt(RMSE_mat, varnames = c("model_scen", "param"), value.name = "relativeRMSE")

df <- df[df$model_scen!="Xiao20_ACi(different Q) & CF",]
df <- df %>%
  separate(model_scen, into = c("model", "datatype"), sep = "_", extra = "merge")

df$model[df$model=="Caemmerer00"] <- "vonC00"
# df$model[df$model=="Busch18"] <- "Busch17"
df$datatype[df$datatype=="ACi(satQ)"] <- "Single A-Ci"
df$datatype[df$datatype=="ACi & AQ"] <- "A-Ci & A-Q"
df$datatype[df$datatype=="ACi(different Q)"] <- "Multiple A-Ci"



newmodels <-  c("FvCB80", "Harley92","vonC00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch18","Xiao20")


datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")

df$param <- factor(df$param,levels=rev(param_name))
df$model <- factor(df$model,levels=newmodels)
df$datatype <- factor(df$datatype,levels=datatypes)


# df$relativeRMSE[df$param=="theta_J" & df$model=="Harley92"] <- NA

midp <- min(df$relativeRMSE,na.rm =T)+(max(df$relativeRMSE,na.rm =T)-min(RMSE_mat,na.rm =T))/2

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

ggplot(df, aes(x = model, y = param, fill = relativeRMSE, label = round(relativeRMSE, 2))) +
  geom_tile(color = "white") +
  # geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
  scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                       limits = c(min(df$relativeRMSE,na.rm =T), max(df$relativeRMSE,na.rm =T)), name = "Relative RMSE (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 11),
        axis.text = element_text(size = 11)) +
  facet_grid(.~datatype) +
  scale_y_discrete(labels = rev(exp_vec))



ggsave(filename = paste0("results/Figures/SFig6_RMSE_all_params.png"),width = 7, height = 7)
