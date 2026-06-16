library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(TeachingDemos)
setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")


load(file="results/Rdata/data_visualization/fig3_MAP_KCO.RData")

load(file="results/Rdata/fig2_sampled_params.RData")
for (i in 1:100){
  sampled_params[[i]]$K_CO <- sampled_params[[i]]$K_C * (1 + 210 / sampled_params[[i]]$K_O)
}

# 1) Extract all true parameter vectors (length 100 each)
true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")
# iter is 1..100, columns are params (Vcmax, Jmax, ...)


param_name <- c("theta_J","alpha_J", "K_CO",
                  "gm","V_tpu","alpha_G", "phi2m", "f_Q","f_pseudo", "h",
                  "g_ch",  "g_wp","max_aG","max_aS","N_max","s","Phi2LL")


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

df$model[df$model=="Caemmerer00"] <- "vonC00"
# df$model[df$model=="Busch18"] <- "Busch17"
df$model[df$model=="Xiao20"] <- "Xiao21"
df$datatype[df$datatype=="ACi(satQ)"] <- "Single A-Ci"
df$datatype[df$datatype=="ACi & AQ"] <- "A-Ci & A-Q"
df$datatype[df$datatype=="ACi(different Q)"] <- "Multiple A-Ci"


spec_cor <- df

spec_cor$param <- factor(spec_cor$param,levels=rev(param_name))

newmodels <-  c("FvCB80", "Harley92","vonC00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch18","Xiao21")
spec_cor$model <- factor(spec_cor$model,levels=newmodels)


datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")

spec_cor$datatype <- factor(spec_cor$datatype,
                            levels=datatypes)

exp_vec <- c(expression(italic(theta[J])),
              expression(italic(alpha[J])),
              expression(italic(K[CO])),
              expression(italic(g[m])),expression(italic(V[TPU])),expression(italic(alpha[G])),
             expression(italic(Phi[IIm])),expression(italic(f[Q])),expression(italic(f[pseudo])),expression(italic("h")),
             expression(italic(g[ch])),expression(italic(g[wp])),
             expression(italic(max[aG])),expression(italic(max[aS])),expression(italic(N[max])),
             expression(italic("s")),expression(italic(Phi[IILL])))

midp <- min(spec_cor$correlation,na.rm =T)+(max(spec_cor$correlation,na.rm =T)-min(cor_mat,na.rm =T))/2



ggplot(spec_cor, aes(x = model, y = param, fill = correlation, label = round(correlation, 2))) +
geom_tile(color = "white") +
geom_text(size = 3, na.rm = TRUE) +  # adds correlation values
scale_fill_gradient2(low = "#77DD77",   mid = "#FFF176",    high = "#FFB347", midpoint = midp,
                     limits = c(min(spec_cor$correlation,na.rm =T), max(spec_cor$correlation,na.rm =T)), name = "Correlation") +
theme_minimal(base_size = 12) +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(size = 12),
      axis.text = element_text(size = 11)) +
facet_grid(.~datatype) + scale_y_discrete(labels = rev(exp_vec))


ggsave(filename = paste0("results/Figures/SFig5_accuracy_remain.png"),width = 8, height = 6)
