library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(ggplot2)

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

generate_curves <- function(model,parameters,model_name,datatype){
  Ci     <- c(50,100,150,200,300,500,800,1000,1200,1500)
  PPFDs  <- c(50, 150, 300, 500, 1200,1800)

  df <- expand.grid(C_i = Ci, PPFD = PPFDs)
  df$O <- 210
  envs <- list(
    C_i  = df$C_i,
    PPFD = df$PPFD,
    O = df$O
  )

  An <- vapply(seq_along(envs$C_i), function(i)
    model(envs$C_i[i], envs$O[i], envs$PPFD[i],
          parameters)[1],
    numeric(1)
  )

  if (model_name == "Xiao20"){

    Phi2 <- vapply(seq_along(envs$C_i), function(i)
      model(envs$C_i[i], envs$O[i], envs$PPFD[i],
            parameters)[2],
      numeric(1)
    )
  }else{Phi2=NA_real_}

  sim_res0 <- list(Ci=df$C_i, PPFD=df$PPFD, An=An, Phi2=Phi2)

  if (datatype=="ACi(satQ)"){ #
    sim_res <- list(Ci=sim_res0$Ci[sim_res0$PPFD==1800],
                    PPFD=sim_res0$PPFD[sim_res0$PPFD==1800],
                    An=sim_res0$An[sim_res0$PPFD==1800])

  }else if(datatype == "ACi & AQ"){#
    sim_res <- list(Ci=c(sim_res0$Ci[sim_res0$PPFD==1800],sim_res0$Ci[sim_res0$Ci==300]),
                    PPFD=c(sim_res0$PPFD[sim_res0$PPFD==1800],sim_res0$PPFD[sim_res0$Ci==300]),
                    An=c(sim_res0$An[sim_res0$PPFD==1800],sim_res0$An[sim_res0$Ci==300]))

  }else if(datatype == "ACi(different Q)"){ #
    sim_res <- list(Ci=sim_res0$Ci,
                    PPFD=sim_res0$PPFD,
                    An=sim_res0$An)
  }else if (datatype == "ACi(different Q) & CF"){ #
    sim_res <- sim_res0
  }

  return(sim_res)
}

mse      <- function(y, yhat) mean((y - yhat)^2, na.rm = TRUE)
rmse     <- function(y, yhat) sqrt(mse(y, yhat))

load(file="results/Rdata/fig2_simulated_curves.RData")

start <- T

for (rep_x in 1:100){
  file <- paste0("results/Rdata/fig2_samples_KCO/fig2_identifiability_KCO_i",rep_x,".RData")
  load(file)

  for (i in 1:length(res_list)){
    temp <- res_list[[i]]
    ind <- which.max(temp$lp__)
    map_vec <- temp[ind,]
    rep <- paste0("rep",rep_x,"_")
    rmind <- which(colnames(temp)%in%c("lp__",".draw"))
    map_vec <- map_vec[,-rmind]
    map_pars <- as.list(map_vec)

    case <- names(res_list)[i]
    check <- strsplit(case,"_")
    datatype <- check[[1]][3]
    model <- check[[1]][2]

    est_curve <- generate_curves(rst_functions[[model]],map_pars,model, datatype)
    true_curve <- simulated_curves[[case]]

    temp <- data.frame(rmse=rmse(true_curve$An, est_curve$An),
                       model=model,
                       datatype=datatype,
                       rep=rep_x)

    if (start){
      score_df <- temp
      start <- F
    }else{
      score_df <- rbind(score_df, temp)
    }
  }
}


score_df <- score_df[score_df$datatype!="ACi(different Q) & CF",]
score_df$model[score_df$model=="Caemmerer00"] <- "vonCaemmerer00"
score_df$model[score_df$model=="Busch18"] <- "Busch17"
score_df$model[score_df$model=="Xiao20"] <- "Xiao21"
score_df$datatype[score_df$datatype=="ACi(satQ)"] <- "Single A-Ci"
score_df$datatype[score_df$datatype=="ACi & AQ"] <- "A-Ci & A-Q"
score_df$datatype[score_df$datatype=="ACi(different Q)"] <- "Multiple A-Ci"

datatypes <- c("Single A-Ci", "A-Ci & A-Q", "Multiple A-Ci")
score_df$datatype <- factor(score_df$datatype,levels= datatypes)

newmodels <-  c("FvCB80", "Harley92","vonCaemmerer00",
                "Ethier04", "Yin04","Dubois07",
                "Tholen12", "Busch17","Xiao21")
score_df$model <- factor(score_df$model,levels= newmodels)


cols <- c("#1B9E77", "#D95F02", "#7570B3")

ggplot(score_df, aes(x = datatype, y = rmse, fill = datatype, colour = datatype)) +
  geom_boxplot( alpha = 0.30) +
  stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)),
               vjust = -0.5, colour = "black", size = 3, fontface = "bold") +

  facet_wrap(.~ model, scales = "fixed") +
  labs(x = "", y = "RMSE", colour = "Data", fill = "Data") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.spacing.x = grid::unit(8, "pt"),
        axis.text.x=element_blank(),
        ) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)


ggsave(filename = paste0("results/Figures/SFig6_RMSE_fitting.png"),width = 7, height = 6)

