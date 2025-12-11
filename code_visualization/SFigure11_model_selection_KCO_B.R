library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")

devtools::load_all()

models <- c("FvCB80", "Harley92", "vonCaemmerer00", "Ethier04", "Yin04",
            "Dubois07", "Tholen12", "Busch17", "Xiao21" )

models0 <- c("FvCB80", "Harley92", "Caemmerer00", "Ethier04", "Yin04",
             "Dubois07", "Tholen12", "Busch18", "Xiao20" )

df_final <- list()
figs <- list()
##################################################
for (scena in 1:2){
  ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
  AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")

  common_plot <- intersect(ACi$Plot, AQ$Plot)
  ACi <- ACi[ACi$Plot%in%common_plot,]
  AQ <- AQ[AQ$Plot%in%common_plot,]


  start <- T
  for (rep_x in 1:length(common_plot)){
    file <- paste0("results/Rdata/fig4_plot_KCO_B/fig4_plot_KCO_i",rep_x,".RData")
    load(file)

    plot_id <- common_plot[rep_x]

    if (scena==1){
      idx <- grep("_ACi$", names(res_list))
    }else if(scena==2){
      idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
    }
    res_list <- res_list[idx]

    for (i in 1:length(res_list)){
      df_pars <- res_list[[i]]
      case <- names(res_list)[i]
      check <- strsplit(case,"_")
      model <- models0[match(check[[1]][2], models0)]

      if (nrow(df_pars)==1 & ncol(df_pars)==10){
        temp <- data.frame( AIC =NA,
                            model=models[match(model, models0)],
                            ID=plot_id)
      }else{
        ind <- which.max(df_pars$lp__)
        map_vec <- df_pars[ind,]
        rep <- paste0("Plot",rep_x,"_")
        loglik <- map_vec$log_lik_sum
        rmind <- which(colnames(df_pars)%in%c("lp__",".draw"))
        map_vec <- map_vec[,-rmind]
        map_pars <- as.list(map_vec)

        AIC_val <- 2 * length(map_pars) - 2 * loglik

        temp <- data.frame(AIC= AIC_val,
                           model=models[match(model, models0)],
                           ID=plot_id)
      }
      if (start){
        score_df <- temp
        start <- F
      }else{
        score_df <- rbind(score_df, temp)
      }
    }
  }

  score_df <- score_df[complete.cases(score_df),]
  score_df1 <- score_df
  score_df1$level <- "Plot"

  score_df1 <- score_df1 %>%
    group_by(ID) %>%
    filter(n_distinct(model) == length(unique(score_df1$model))) %>%
    ungroup()

  #####################################################################
  ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
  AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
  ACi$ID <- paste0(ACi$Plot,"_",ACi$Repeat)
  AQ$ID <- paste0(AQ$Plot,"_",ACi$Repeat)

  plant_id <- intersect(ACi$ID, AQ$ID)

  start <- T
  for (rep_x in 1:length(plant_id)){
    file <- paste0("results/Rdata/fig4_plant_KCO_B/fig4_plant_fit_i",rep_x,".RData")
    load(file)

    if (scena==1){
      idx <- grep("_ACi$", names(res_list))
    }else if(scena==2){
      idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
    }
    res_list <- res_list[idx]
    res_models <- sub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", names(res_list))

    plant <- plant_id[rep_x]
    for (i in 1:length(res_list)){
      df_pars <- res_list[[i]]
      model <- models0[match(res_models[i], models0)]
      if (nrow(df_pars)==1 & ncol(df_pars)==10){
        temp <- data.frame( AIC =NA,
                            model=models[match(model, models0)],
                            ID=plant)
      }else{
        ind <- which.max(df_pars$lp__)
        map_vec <- df_pars[ind,]
        loglik <- map_vec$log_lik_sum
        rep <- paste0("Plot",rep_x,"_")
        rmind <- which(colnames(df_pars)%in%c("lp__",".draw"))
        map_vec <- map_vec[,-rmind]
        map_pars <- as.list(map_vec)

        AIC_val <- 2 * length(map_pars) - 2 * loglik

        temp <- data.frame(AIC= AIC_val,
                           model=models[match(model, models0)],
                           ID=plant)
      }

      if (start){
        score_df <- temp
        start <- F
      }else{
        score_df <- rbind(score_df, temp)
      }
    }
  }
  score_df <- score_df[complete.cases(score_df),]
  score_df2 <- score_df
  score_df2$level <- "Plant"
  score_df2 <- score_df2 %>%
    group_by(ID) %>%
    filter(n_distinct(model) == length(unique(score_df2$model))) %>%
    ungroup()
  #####################################################################
  ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
  AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
  ref <- read.csv("data/2022_barley_reference.csv")

  AQ <- AQ[AQ$Flag_removal!="x",]
  ACi$genotype <- ref$Accession[match(ACi$Plot,ref$PlotID)]
  AQ$genotype <- ref$Accession[match(AQ$Plot,ref$PlotID)]

  genotype_id <- intersect(ACi$genotype, AQ$genotype)

  start <- T
  for (rep_x in 1:length(genotype_id)){
    file <- paste0("results/Rdata/fig4_genotype_KCO_B/fig4_genotype_fit_i",rep_x,".RData")
    load(file)
    geno <- genotype_id[rep_x]
    if (scena==1){
      idx <- grep("_ACi$", names(res_list))
    }else if(scena==2){
      idx <- grep("ACi.*AQ|AQ.*ACi", names(res_list))
    }
    res_list <- res_list[idx]
    res_models <- sub("^[^_]+_([^_]+)_.*$", "\\1", names(res_list))

    for (i in 1:length(res_list)){
      df_pars <- res_list[[i]]
      model <- models0[match(res_models[i], models0)]
      if (nrow(df_pars)==1 & ncol(df_pars)==10){
        temp <- data.frame( AIC =NA,
                            model=models[match(model, models0)],
                            ID=geno)
      }else{
        ind <- which.max(df_pars$lp__)
        map_vec <- df_pars[ind,]
        loglik <- map_vec$log_lik_sum
        rmind <- which(colnames(df_pars)%in%c("lp__",".draw"))
        map_vec <- map_vec[,-rmind]
        map_pars <- as.list(map_vec)

        AIC_val <- 2 * length(map_pars) - 2 * loglik

        temp <- data.frame(AIC= AIC_val,
                           model=models[match(model, models0)],
                           ID=geno)
      }

      if (start){
        score_df <- temp
        start <- F
      }else{
        score_df <- rbind(score_df, temp)
      }
    }
  }
  score_df <- score_df[complete.cases(score_df),]
  score_df3 <- score_df
  score_df3$level <- "Genotype"

  score_df3 <- score_df3 %>%
    group_by(ID) %>%
    filter(n_distinct(model) == length(unique(score_df3$model))) %>%
    ungroup()
  #####################################################################
  score_df <- rbind(score_df1,score_df2)
  score_df <- rbind(score_df,score_df3)
  score_df$model <- factor(score_df$model,
                           levels= models)

  check <- score_df[is.na(score_df$AIC),]
  toremove <- unique(check$ID)
  ind <- which(score_df$ID%in%toremove)
  score_df <- score_df[setdiff(1:nrow(score_df),ind),]


  score_df$level <- factor(score_df$level, levels=c("Plant","Plot", "Genotype"))

  df_final[[scena]] <- score_df

  figs[[scena]] <- ggplot(score_df, aes(x = level, y = AIC, fill = level, colour = level)) +
    geom_boxplot( alpha = 0.30) +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)),
                 vjust = -0.5, colour = "black", size = 3.5) + #, fontface = "bold"

    facet_wrap(.~ model, scales = "fixed") +
    labs(x = "", y = "AIC", colour = "Data", fill = "Data") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.spacing.x = grid::unit(8, "pt"),
          axis.text.x=element_blank()) +
    coord_cartesian(ylim=c(0, 400))
}


save(df_final, file="results/Rdata/data_visualization/SFig11.RData")

########################################################################
library("ggpubr")
legs <- c("none","bottom")
for (scena in 1:2){
figs[[scena]] <- ggplot(df_final[[scena]], aes(x = level, y = AIC, fill = level, colour = level)) +
  geom_boxplot( alpha = 0.30) +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 1)),
               vjust = -0.5, colour = "black", size = 3.5) + #, fontface = "bold"

  facet_wrap(.~ model, scales = "fixed") +
  labs(x = "", y = "AIC", colour = "", fill = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = legs[scena],
        legend.text = element_text(size=11),
        panel.spacing.x = grid::unit(8, "pt"),
        axis.text.x=element_blank(),
        strip.text = element_text(size = 10)) +
  coord_cartesian(ylim=c(0, 700))
}
ggarrange(figs[[1]],figs[[2]],labels = c('a', 'b'), ncol=1, heights = c(0.8,1))


ggsave(filename = paste0("results/Figures/SFig11_model_AIC.png"),width = 7, height = 7)

