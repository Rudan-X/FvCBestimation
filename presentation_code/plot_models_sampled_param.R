library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(reshape2)

devtools::load_all()
# check check
# Grid search for environmental input:
# Ci <- c(20,50,100,150,200,300,400,500,600,700,800,1000,1200)
# PPFDs <- c(20,100,300,500,700,1000,1200,1500)

Ci     <- c(50,100,150,200,300,500,600,700,800,1000,1200)
PPFDs  <- c(1800) # 50,150,300,500,1100,

sampling_res <- readRDS("results/sampling_scenario1_100000.rds")


model_names <- c("FvCB (1980) (Wc,Wj) (Ci~=Cc)", "Harley (1992) (Ac,Aj) (Ci=Cc-A/gm)",
                 "von Caemmerer (2000) (Wc,Wj,Wp) (Ci~=Cc)","Ethier & Long (2004) (Ac,Aj) (Ci=Cc-A/gm)",
                 "Yin (2004) (Wc,Wj) (Ci~=Cc)","Dubois (2007) (Ac,Aj,Ap) (Ci~=Cc)",
                 "Tholen (2012) (Wc,Wj) (Ci=Cc-A/gm)","Busch (2018) (Wc,Wj,Wp) (Ci~=Cc)")

for (m in 1:length(model_name)){
  temp <- sampling_res$simA[[m]][,6,]
  rem.ind <- union(which(temp[2,] < (-2)),which(temp[3,] < (-2)))
  rem.ind <- union(rem.ind, which(temp[4,] < (-2)))
  df0 <- melt(temp[, -rem.ind])

  df0$Model <- model_names[m]
  if (m==1){
    df_final <- df0
  }else{
    df_final <- rbind(df_final,df0)
  }
}

summary_df <- df_final %>%
  group_by(Model, Ci) %>%
  summarise(
    Median_A =round(median(value, na.rm = TRUE),2),
    Std_A = round(sd(value, na.rm = TRUE),3),
    lower_A = quantile(value, 0.25, na.rm = TRUE),
    upper_A = quantile(value, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

new_name <- c("FvCB (1980) (Wc,Wj) (Ci~=Cc)", "Yin (2004) (Wc,Wj) (Ci~=Cc)",
                "Harley (1992) (Ac,Aj) (Ci=Cc-A/gm)", "Ethier & Long (2004) (Ac,Aj) (Ci=Cc-A/gm)",
                "Tholen (2012) (Wc,Wj) (Ci=Cc-A/gm)",
                "von Caemmerer (2000) (Wc,Wj,Wp) (Ci~=Cc)","Dubois (2007) (Ac,Aj,Ap) (Ci~=Cc)",
                "Busch (2018) (Wc,Wj,Wp) (Ci~=Cc)")

summary_df$Model <- factor(summary_df$Model, levels = new_name)


ggplot(summary_df, aes(x = Ci)) +

  geom_errorbar(aes(ymin = lower_A, ymax = upper_A), width = 1, color = "orange") +
  geom_point(aes(y = Median_A),shape = 1, size = 2.5, color = "purple") +

  facet_wrap(~Model,ncol=4,dir = "v") +
  labs(
    x = expression(C[i]~"(µmol mol"^{-1}*")"),
    y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
    color = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-5,60))


ggsave("results/Figures/Sampled_parameter.png",width = 12, height = 7)


#
#
# df_final$Model <- factor(df_final$Model, levels = model_name)
#
# ggplot(df_final, aes(x = Ci, y = value)) +
#   geom_point(outlier.shape = 1) +
#   facet_wrap(~Model,ncol=2)+
#   coord_cartesian(ylim = c(-5,60))+
# theme_bw() +
#   ylab("A") +
#   xlab("Ci")
#
#
# ggplot(df_final, aes(x = Ci, y = value, group = sample)) +
#   geom_line(alpha = 0.3) +
#   facet_wrap(~Model,ncol=2) +
#   labs(
#     x = expression(C[i]~"(µmol mol"^{-1}*")"),
#     y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
#     color = "Model"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim = c(-5,60))
#
