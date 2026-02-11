library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)



filen<-paste0("R/")
myls <- list.files(path=filen,pattern="*.R")
myls<-paste0(filen,myls)
invisible(sapply(myls,FUN=source))



devtools::load_all()


load(file="results/Rdata/fig2_simulated_curves.RData")

rst_FvCB <- stan_model(model_code = stan_FvCB80_B)
rst_Harley <- stan_model(model_code = stan_Harley92_B)
rst_Caemmerer <- stan_model(model_code = stan_Caemmerer00_B)
rst_Ethier <- stan_model(model_code = stan_Ethier04_B)
rst_Yin <- stan_model(model_code = stan_Yin04_B)
rst_Dubois <- stan_model(model_code = stan_Dubois07_B)
rst_Tholen <- stan_model(model_code = stan_Tholen12_B)
rst_Busch <- stan_model(model_code = stan_Busch18_B)
rst_Xiao_A <- stan_model(model_code = stan_Xiao20_An_B)
rst_Xiao_Phi <- stan_model(model_code = stan_Xiao20_Phi_B)

rst_models <- list(
  FvCB80 = rst_FvCB,
  Harley92 = rst_Harley,
  Caemmerer00  = rst_Caemmerer,
  Ethier04  = rst_Ethier,
  Yin04 = rst_Yin,
  Dubois07 = rst_Dubois,
  Tholen12 = rst_Tholen,
  Busch18 = rst_Busch,
  Xiao20 = list(rst_Xiao_A,rst_Xiao_Phi)
)


param_candidates <- c("V_cmax","J_max","R_d","gamma_star",
                      "theta_J","alpha_J", "K_C", "K_O","V_tpu", "alpha_G","gm",
                      "phi2m", "f_Q", "f_pseudo", "h", "N_max", "max_aG", "max_aS",
                      "g_ch","g_wp","s","Phi2LL")

load(file="results/Rdata/fig2_sampled_params.RData")

true_params_df <- map_dfr(sampled_params, ~as.data.frame(.x), .id = "iter")[,-1]

median_params <- apply(true_params_df, 2, median, na.rm = TRUE)
median_params <- as.list(median_params)

build_stan_data <- function(sim_res) {
  sim_res <- sim_res[complete.cases(sim_res),]
  true_val <- median_params
  keep <- setdiff(names(true_val),c("V_cmax","J_max","R_d","gamma_star"))
  sampled_val <- true_val[keep]
  n=length(sim_res$Ci)
  data <- list(
    N     = as.integer(n),
    Ci    = as.numeric(sim_res$Ci),
    O     = as.numeric(sim_res$O),
    PPFD  = as.numeric(sim_res$PPFD),
    A_obs = as.numeric(sim_res$An),
    has_Cc     = 0L,
    Cc         = rep(0.0, n))

  data <- c(data,sampled_val)
  if ("Phi2"%in%colnames(sim_res)){
    data$Phi2_obs <- as.numeric(sim_res$Phi2)
  }
  return(data)
}


datatypes <- c("ACi(satQ)", "ACi & AQ", "ACi(different Q)", "ACi(different Q) & CF")


fit_one_repitition <- function(rep_id,rst_models,
                               iter = 3000, warmup = 1000, chains = 4,
                               adapt_delta = 0.95, max_treedepth = 12) {

  res_list <- list()
  # storage
  for (m in 1:length(rst_models)) {

    if (m==9){
      dl <- 4
    }else{dl <- 3}

    for (d in 1:dl){
      if (m==9){
        if (d!=4){
          stan_m <- rst_models[[m]][[1]]
        }else{
          stan_m <- rst_models[[m]][[2]]
        }
      }else{
        stan_m <- rst_models[[m]]
      }

      print(paste0("rep",rep_id,"_",names(rst_models)[m],"_",datatypes[d]))

      case_name <- paste0("rep",rep_id,"_",names(rst_models)[m],"_",datatypes[d])
      x <- which(names(simulated_curves)==case_name)
      sim_res <- simulated_curves[[x]]
      sim_res$O <- 210

      sdata <- build_stan_data(sim_res)
      # fit_map <- optimizing(stan_m, data = sdata)
      fit <- sampling(
        stan_m, data = sdata, cores=4,
        chains = chains, iter = iter, warmup = warmup,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        refresh = 0
      )

      dr <- posterior::as_draws_df(fit)
      keep <- intersect(colnames(dr), param_candidates)

      if (m==9){
        if (d!=4){
          keep <- c(keep, "sigma_A","lp__",".draw")
        }else{
          keep <- c(keep, "sigma_A", "sigma_Phi","lp__",".draw")
        }
      }else{
        keep <- c(keep, "sigma","lp__",".draw")
      }


      res_list[[length(res_list) + 1]] <- dr[, keep, drop = FALSE]
      names(res_list)[length(res_list)] <- paste0("rep",rep_x,"_",names(rst_models)[m],"_",datatypes[d])
    }
  }
  save(res_list, file=paste0("results/Rdata/fig3_samples_median_fixed/fig3_core_estimate_median_fixed_i",rep_id,".RData"))
}






for (rep_x in 76:100){
  fit_one_repitition(rep_x, rst_models = rst_models)
}


for (rep_x in 1:25){
  fit_one_repitition(rep_x, rst_models = rst_models)
}
for (rep_x in 26:50){
  fit_one_repitition(rep_x, rst_models = rst_models)
}


for (rep_x in 51:75){
  fit_one_repitition(rep_x, rst_models = rst_models)
}
