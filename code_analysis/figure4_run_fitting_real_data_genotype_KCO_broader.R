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

ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
ref <- read.csv("data/2022_barley_reference.csv")

AQ <- AQ[AQ$Flag_removal!="x",]
ACi$genotype <- ref$Accession[match(ACi$Plot,ref$PlotID)]
AQ$genotype <- ref$Accession[match(AQ$Plot,ref$PlotID)]

genotype_id <- intersect(ACi$genotype, AQ$genotype)

ACi <- ACi[ACi$genotype%in%genotype_id,]
AQ <- AQ[AQ$genotype%in%genotype_id,]


rst_FvCB <- stan_model(model_code = stan_FvCB80_D)
rst_Harley <- stan_model(model_code = stan_Harley92_D)
rst_Caemmerer <- stan_model(model_code = stan_Caemmerer00_D)
rst_Ethier <- stan_model(model_code = stan_Ethier04_D)
rst_Yin <- stan_model(model_code = stan_Yin04_D)
rst_Dubois <- stan_model(model_code = stan_Dubois07_D)
rst_Tholen <- stan_model(model_code = stan_Tholen12_D)
rst_Busch <- stan_model(model_code = stan_Busch18_D)
rst_Xiao_A <- stan_model(model_code = stan_Xiao20_An_D)



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

####

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



fit_one_repitition <- function(ind_x,rst_models,   # list of compiled Stan models
                               iter = 3000, warmup = 1000, chains = 4,
                               adapt_delta = 0.95, max_treedepth = 12) {
  datatypes <- c("ACi","ACi & AQ")
  param_candidates <- c("V_cmax","J_max","R_d","gamma_star",
                        "theta_J","alpha_J", "K_CO","V_tpu", "alpha_G","gm",
                        "phi2m", "f_Q", "f_pseudo", "h", "N_max", "max_aG", "max_aS",
                        "g_ch","g_wp","s","Phi2LL")
  res_list <- list()
  p_id <- genotype_id[ind_x]
  for (m in 1:length(rst_models)) {
    stan_m <- rst_models[[m]]
    for (d in 1:2){
      if (d==1){
        curve <- ACi[ACi$genotype==p_id,c("Ci","PPFD","Photo")]
        colnames(curve)[3] <- "An"
      }else{
        curve1 <- ACi[ACi$genotype==p_id,c("Ci","PPFD","Photo")]
        colnames(curve1)[3] <- "An"
        curve2 <- AQ[AQ$genotype==p_id,c("Ci","PAR","A")]
        colnames(curve2)[c(2,3)] <- c("PPFD","An")
        curve <- rbind(curve1,curve2)
      }
      curve$O <- 210

      print(paste0("Genotype",ind_x,"_",names(rst_models)[m],"_",datatypes[d]))

      case_name <- paste0("Genotype",p_id,"_",names(rst_models)[m],"_",datatypes[d])

      sdata <- build_stan_data(curve)

      fit <- try(
        sampling(
          stan_m, data = sdata,
          chains = chains, iter = iter, warmup = warmup, cores = 4,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
          refresh = 0
        ),
        silent = TRUE
      )

      if (length(rstan::extract(fit, permuted = FALSE)) == 0) {
        message("Model fitting failed for this case.")
        fit <- NULL
        res_list[[length(res_list) + 1]] <- matrix(NA, 1,10)
      }else{
        dr <- posterior::as_draws_df(fit)
        keep <- intersect(colnames(dr), param_candidates)
        keep <- c(keep, "sigma", "log_lik_sum", "lp__",".draw")
        dr <- dr[, keep, drop = FALSE]
        res_list[[length(res_list) + 1]] <- dr

      }
      names(res_list)[length(res_list)] <- case_name
    }
  }
  save(res_list, file=paste0("results/Rdata/fig4_genotype_KCO_B/fig4_genotype_fit_i",ind_x,".RData"))
}

for (rep_x in 41:85){
  fit_one_repitition(rep_x, rst_models = rst_models)
}



for (rep_x in 1:40){
  fit_one_repitition(rep_x, rst_models = rst_models)
}




