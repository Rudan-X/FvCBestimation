library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")

filen<-paste0("R/")
myls <- list.files(path=filen,pattern="*.R")
myls<-paste0(filen,myls)
invisible(sapply(myls,FUN=source))



devtools::load_all()

# load(file="results/Rdata/fig2_sampled_params.RData")

load(file="results/Rdata/fig2_simulated_curves.RData")

rst_FvCB <- stan_model(model_code = stan_FvCB80_C)
rst_Harley <- stan_model(model_code = stan_Harley92_C)
rst_Caemmerer <- stan_model(model_code = stan_Caemmerer00_C)
rst_Ethier <- stan_model(model_code = stan_Ethier04_C)
rst_Yin <- stan_model(model_code = stan_Yin04_C)
rst_Dubois <- stan_model(model_code = stan_Dubois07_C)
rst_Tholen <- stan_model(model_code = stan_Tholen12_C)
rst_Busch <- stan_model(model_code = stan_Busch18_C)
rst_Xiao_A <- stan_model(model_code = stan_Xiao20_An_C)


rst_models <- list(
  FvCB80 = rst_FvCB,
  Harley92 = rst_Harley,
  Caemmerer00  = rst_Caemmerer,
  Ethier04  = rst_Ethier,
  Yin04 = rst_Yin,
  Dubois07 = rst_Dubois,
  Tholen12 = rst_Tholen,
  Busch18 = rst_Busch,
  Xiao20 = rst_Xiao_A)


param_candidates <- c("V_cmax","J_max","R_d","gamma_star",
                      "theta_J","alpha_J", "K_CO","V_tpu", "alpha_G","gm",
                      "phi2m", "f_Q", "f_pseudo", "h", "N_max", "max_aG", "max_aS",
                      "g_ch","g_wp","s","Phi2LL")


build_stan_data <- function(sim_res) {
  sim_res <- sim_res[complete.cases(sim_res),]
  n=length(sim_res$Ci)
  data <- list(
    N     = as.integer(n),
    Ci    = as.numeric(sim_res$Ci),
    O     = as.numeric(sim_res$O),
    PPFD  = as.numeric(sim_res$PPFD),
    A_obs = as.numeric(sim_res$An),
    has_Cc     = 0L,
    Cc         = rep(0.0, n)
  )
  if ("Phi2"%in%colnames(sim_res)){
    data$Phi2_obs <- as.numeric(sim_res$Phi2)
  }
  return(data)
}



datatypes <- c("ACi(satQ)", "ACi & AQ", "ACi(different Q)") #, "ACi(different Q) & CF"


fit_one_repitition <- function(rep_id,rst_models,   # list of compiled Stan models
                               iter = 3000, warmup = 1000, chains = 4,
                               adapt_delta = 0.95, max_treedepth = 12) {

  res_list <- list()
  # storage
  for (m in 1:length(rst_models)) {

    for (d in 1:3){

      stan_m <- rst_models[[m]]

      print(paste0("rep",rep_x,"_",names(rst_models)[m],"_",datatypes[d]))

      case_name <- paste0("rep",rep_id,"_",names(rst_models)[m],"_",datatypes[d])
      x <- which(names(simulated_curves)==case_name)
      sim_res <- simulated_curves[[x]]
      sim_res$O <- 210

      sdata <- build_stan_data(sim_res)
      fit <- sampling(
        stan_m, data = sdata,
        chains = chains, cores=4,iter = iter, warmup = warmup,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        refresh = 0
      )
      dr <- posterior::as_draws_df(fit)
      keep <- intersect(colnames(dr), param_candidates)
      keep <- c(keep, "lp__",".draw")
      dr <- dr[, keep, drop = FALSE]

      res_list[[length(res_list) + 1]] <- dr
      names(res_list)[length(res_list)] <- paste0("rep",rep_x,"_",names(rst_models)[m],"_",datatypes[d])
    }
  }
  save(res_list, file=paste0("results/Rdata/fig2_samples_KCO/fig2_identifiability_KCO_i",rep_id,".RData"))
}

for (rep_x in 77:100){
  fit_one_repitition(rep_x, rst_models = rst_models)
}

for (rep_x in 27:50){
  fit_one_repitition(rep_x, rst_models = rst_models)
}


for (rep_x in 2:26){
  fit_one_repitition(rep_x, rst_models = rst_models)
}





for (rep_x in 52:76){
  fit_one_repitition(rep_x, rst_models = rst_models)
}



