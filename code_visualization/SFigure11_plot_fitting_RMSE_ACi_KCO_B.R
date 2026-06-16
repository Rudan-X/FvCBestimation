library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")


devtools::load_all()

rst_functions <- stan_model(model_code = stan_all_functions_KCO)
expose_stan_functions(rst_functions)


call_function_FvCB80 <- function(Ci, O, PPFD, pars) {
  function_FvCB80(Ci, 0L, 0, O, PPFD,
                  pars$V_cmax, pars$J_max, pars$R_d,
                  pars$K_CO, pars$gamma_star,
                  pars$alpha_J, pars$theta_J)
}

call_function_Harley92 <- function(Ci, O, PPFD, pars) {
  function_Harley92(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_CO, pars$gamma_star,
                    pars$alpha_J,pars$gm)
}

call_function_Caemmerer00 <- function(Ci, O, PPFD, pars) {
  function_Caemmerer00(Ci, 0L, 0, O, PPFD,
                       pars$V_cmax, pars$J_max, pars$R_d,
                       pars$K_CO, pars$gamma_star,
                       pars$alpha_J, pars$theta_J,
                       pars$V_tpu, pars$alpha_G)
}

call_function_Ethier04 <- function(Ci, O, PPFD, pars) {
  function_Ethier04(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_CO, pars$gamma_star,
                    pars$alpha_J, pars$theta_J,
                    pars$gm)
}

call_function_Yin04 <- function(Ci, O, PPFD, pars) {
  function_Yin04(Ci, 0L, 0, O, PPFD,
                 pars$V_cmax, pars$J_max, pars$R_d,
                 pars$K_CO, pars$gamma_star,
                 pars$phi2m, pars$theta_J,
                 pars$f_Q, pars$f_pseudo, pars$h)
}


call_function_Dubois07 <- function(Ci, O, PPFD, pars) {
  function_Dubois07(Ci, 0L, 0, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_CO, pars$gamma_star,
                    pars$alpha_J, pars$theta_J,
                    pars$V_tpu, pars$alpha_G)
}


call_function_Tholen12 <- function(Ci, O, PPFD, pars) {
  function_Tholen12(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_CO, pars$gamma_star,
                    pars$alpha_J,pars$theta_J, pars$g_ch ,pars$g_wp)
}


call_function_Busch18 <- function(Ci, O, PPFD, pars) {
  function_Busch18(Ci, 0L, 0, O, PPFD,
                   pars$V_cmax, pars$J_max, pars$R_d,
                   pars$K_CO, pars$gamma_star,
                   pars$alpha_J, pars$theta_J,
                   pars$V_tpu, pars$N_max, pars$max_aG, pars$max_aS)
}

call_function_Xiao20 <- function(Ci, O, PPFD, pars) {
  function_Xiao20(Ci, O, PPFD,
                  pars$V_cmax, pars$J_max, pars$R_d,
                  pars$K_CO, pars$gamma_star,
                  pars$s, pars$Phi2LL, pars$theta_J,
                  pars$gm)
}

rst_functions <- list(
  FvCB80     = call_function_FvCB80,
  Harley92 = call_function_Harley92,
  Caemmerer00  = call_function_Caemmerer00,
  Ethier04   = call_function_Ethier04,
  Yin04 = call_function_Yin04,
  Dubois07   = call_function_Dubois07,
  Tholen12 = call_function_Tholen12,
  Busch18 = call_function_Busch18,
  Xiao20 = call_function_Xiao20
)

generate_curves <- function(model,parameters,Ci, PPFDs){
  envs <- list(
    C_i  = Ci,
    PPFD = PPFDs,
    O = rep(210,length(Ci))
  )

  An <- vapply(seq_along(envs$C_i), function(i)
    model(envs$C_i[i], envs$O[i], envs$PPFD[i],
          parameters)[1],
    numeric(1)
  )
  return(An)
}

mse      <- function(y, yhat) mean((y - yhat)^2, na.rm = TRUE)
rmse     <- function(y, yhat) sqrt(mse(y, yhat))


models <- c("FvCB80", "Harley92", "vonC00", "Ethier04", "Yin04",
            "Dubois07", "Tholen12", "Busch18", "Xiao21" )

models0 <- c("FvCB80", "Harley92", "Caemmerer00", "Ethier04", "Yin04",
             "Dubois07", "Tholen12", "Busch18", "Xiao20" )

df_final <- list()
setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation_not_uploaded_true/")
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


    curve1 <- ACi[ACi$Plot==plot_id,c("Ci","PPFD","Photo")]
    colnames(curve1)[3] <- "An"
    curve2 <- AQ[AQ$Plot==plot_id,c("Ci","PAR","A")]
    colnames(curve2)[c(2,3)] <- c("PPFD","An")
    if (scena==1){
      curve <- curve1
    }else{
      curve <- rbind(curve1,curve2)
    }

    for (i in 1:length(res_list)){
      df_pars <- res_list[[i]]
      case <- names(res_list)[i]
      check <- strsplit(case,"_")
      model <- models0[match(check[[1]][2], models0)]

      if (nrow(df_pars)==1 & ncol(df_pars)==10){
        temp <- data.frame( rmse =NA,
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

        est_curve <- try(
          generate_curves(rst_functions[[model]], map_pars, curve$Ci, curve$PPFD),
          silent = TRUE
        )

        if (inherits(est_curve, "try-error")){est_curve <- NA}
        temp <- data.frame(rmse=rmse(curve$An,est_curve),
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

    curve1 <- ACi[ACi$ID==plant,c("Ci","PPFD","Photo")]
    colnames(curve1)[3] <- "An"
    curve2 <- AQ[AQ$ID==plant,c("Ci","PAR","A")]
    colnames(curve2)[c(2,3)] <- c("PPFD","An")

    if (scena==1){
      curve <- curve1
    }else{
      curve <- rbind(curve1,curve2)
    }

    for (i in 1:length(res_list)){
      df_pars <- res_list[[i]]
      model <- models0[match(res_models[i], models0)]
      if (nrow(df_pars)==1 & ncol(df_pars)==10){
        temp <- data.frame( rmse =NA,
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

        est_curve <- try(
          generate_curves(rst_functions[[model]], map_pars, curve$Ci, curve$PPFD),
          silent = TRUE
        )
        if (inherits(est_curve, "try-error")){est_curve <- NA}

        temp <- data.frame(rmse=rmse(curve$An,est_curve),
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

    curve1 <- ACi[ACi$genotype==geno,c("Ci","PPFD","Photo")]
    colnames(curve1)[3] <- "An"
    curve2 <- AQ[AQ$genotype==geno,c("Ci","PAR","A")]
    colnames(curve2)[c(2,3)] <- c("PPFD","An")
    if (scena==1){
      curve <- curve1
    }else{
      curve <- rbind(curve1,curve2)
    }

    for (i in 1:length(res_list)){
      df_pars <- res_list[[i]]
      model <- models0[match(res_models[i], models0)]
      if (nrow(df_pars)==1 & ncol(df_pars)==10){
        temp <- data.frame( rmse =NA,
                            model=models[match(model, models0)],
                            ID=geno)
      }else{
        ind <- which.max(df_pars$lp__)
        map_vec <- df_pars[ind,]
        loglik <- map_vec$log_lik_sum
        rmind <- which(colnames(df_pars)%in%c("lp__",".draw"))
        map_vec <- map_vec[,-rmind]
        map_pars <- as.list(map_vec)

        est_curve <- try(
          generate_curves(rst_functions[[model]], map_pars, curve$Ci, curve$PPFD),
          silent = TRUE
        )
        if (inherits(est_curve, "try-error")){est_curve <- NA}

        temp <- data.frame(rmse=rmse(curve$An,est_curve),
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

  check <- score_df[is.na(score_df$rmse),]
  toremove <- unique(check$ID)
  ind <- which(score_df$ID%in%toremove)
  score_df <- score_df[setdiff(1:nrow(score_df),ind),]


  score_df$level <- factor(score_df$level, levels=c("Plant","Plot", "Genotype"))

  df_final[[scena]] <- score_df
}

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")
save(df_final, file="results/Rdata/data_visualization/SFig10_RMSE.RData")

load(file="results/Rdata/data_visualization/SFig10_RMSE.RData")


library("ggpubr")
figs <- list()

legs <- c("none","bottom")
for (scena in 1:2){
  figs[[scena]] <- ggplot(df_final[[scena]], aes(x = level, y = rmse, fill = level, colour = level)) +
    geom_boxplot( alpha = 0.30) +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)),
                 vjust = -1, colour = "black", size = 3.5) + #, fontface = "bold"

    facet_wrap(.~ model, scales = "fixed") +
    labs(x = "", y = "RMSE", colour = "", fill = "") +
    theme_minimal(base_size = 12) +
    theme(legend.position = legs[scena],
          legend.text = element_text(size=11),
          panel.spacing.x = grid::unit(8, "pt"),
          axis.text.x=element_blank(),
          strip.text = element_text(size = 10))+
    coord_cartesian(ylim=c(0, 15))
}
ggarrange(figs[[1]],figs[[2]],labels = c('a', 'b'), ncol=1, heights = c(0.8,1))


ggsave(filename = paste0("results/Figures/SFig11_model_RMSE.png"),width = 7, height = 7)

