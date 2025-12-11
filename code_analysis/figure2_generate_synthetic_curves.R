library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)


devtools::load_all()

sample_true_pars <- function(){
  list(
    V_cmax     = rtn(1, 50, 200),
    J_max      = rtn(1, 100, 250),
    R_d        = rtn(1, 0.5, 5),
    gamma_star = rtn(1, 10, 35),
    V_tpu      = rtn(1, 8, 20),
    gm         = rtn(1, 0.1, 0.7),
    alpha_J    = rtn(1, 0.3, 0.5),
    theta_J    = rtn(1, 0.5, 0.9),
    K_C        = rtn(1, 200, 300),
    K_O        = rtn(1, 150, 300),
    alpha_G    = rtn(1, 0.3, 0.8),
    phi2m      = rtn(1, 0.5, 0.9),
    f_Q        = rtn(1, 0.5, 1),
    f_pseudo   = rtn(1, 0, 0.3),
    h          = rtn(1, 3, 14/3),
    max_aG     = rtn(1, 0.01, 0.2),
    max_aS     = rtn(1, 0.01, 0.5),
    N_max      = rtn(1, 0.1, 2),
    g_ch       = rtn(1, 0.05, 0.7),
    g_wp       = rtn(1, 0.1, 2),
    s          = rtn(1, 0.3, 0.6)
  )
}
# sampled_params <- list()
# for (i in 1:100){
#   sampled_params[[i]] <- sample_true_pars()
# }
#
# save(sampled_params, file="results/Rdata/fig2_sampled_params.RData")
#

load(file="results/Rdata/fig2_sampled_params.RData")

# for (i in 1:100){
#   temp <- sample_true_pars()
#   while (temp$alpha_J/temp$s>1){
#     temp <- sample_true_pars()
#   }
#   sampled_params[[i]]$s <-temp$s
#   sampled_params[[i]]$Phi2LL <- temp$alpha_J/temp$s
# }
#
# save(sampled_params, file="results/Rdata/fig2_sampled_params.RData")

rst_functions <- stan_model(model_code = stan_all_functions)
expose_stan_functions(rst_functions)


call_function_FvCB80 <- function(Ci, O, PPFD, pars) {
  function_FvCB80(Ci, 0L, 0, O, PPFD,
                  pars$V_cmax, pars$J_max, pars$R_d,
                  pars$K_C, pars$K_O, pars$gamma_star,
                  pars$alpha_J, pars$theta_J)
}

call_function_Harley92 <- function(Ci, O, PPFD, pars) {
  function_Harley92(Ci, O, PPFD,
                  pars$V_cmax, pars$J_max, pars$R_d,
                  pars$K_C, pars$K_O, pars$gamma_star,
                  pars$alpha_J,pars$gm)
}

call_function_Caemmerer00 <- function(Ci, O, PPFD, pars) {
  function_Caemmerer00(Ci, 0L, 0, O, PPFD,
                       pars$V_cmax, pars$J_max, pars$R_d,
                       pars$K_C, pars$K_O, pars$gamma_star,
                       pars$alpha_J, pars$theta_J,
                       pars$V_tpu, pars$alpha_G)
}

call_function_Ethier04 <- function(Ci, O, PPFD, pars) {
  function_Ethier04(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
                    pars$alpha_J, pars$theta_J,
                    pars$gm)
}

call_function_Yin04 <- function(Ci, O, PPFD, pars) {
  function_Yin04(Ci, 0L, 0, O, PPFD,
                 pars$V_cmax, pars$J_max, pars$R_d,
                 pars$K_C, pars$K_O, pars$gamma_star,
                 pars$phi2m, pars$theta_J,
                 pars$f_Q, pars$f_pseudo, pars$h)
}


call_function_Dubois07 <- function(Ci, O, PPFD, pars) {
  function_Dubois07(Ci, 0L, 0, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
                    pars$alpha_J, pars$theta_J,
                    pars$V_tpu, pars$alpha_G)
}


call_function_Tholen12 <- function(Ci, O, PPFD, pars) {
  function_Tholen12(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
                    pars$alpha_J,pars$theta_J, pars$g_ch ,pars$g_wp)
}


call_function_Busch18 <- function(Ci, O, PPFD, pars) {
  function_Busch18(Ci, 0L, 0, O, PPFD,
                   pars$V_cmax, pars$J_max, pars$R_d,
                   pars$K_C, pars$K_O, pars$gamma_star,
                   pars$alpha_J, pars$theta_J,
                   pars$V_tpu, pars$N_max, pars$max_aG, pars$max_aS)
}

call_function_Xiao20 <- function(Ci, O, PPFD, pars) {
  function_Xiao20(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
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




param_candidates <- c("V_cmax","J_max","R_d","gamma_star",
                      "theta_J","alpha_J", "K_C", "K_O","V_tpu", "alpha_G","gm",
                      "phi2m", "f_Q", "f_pseudo", "h", "N_max", "max_aG", "max_aS",
                      "s","Phi2LL")

generate_curves <- function(model,parameters,model_name){
  Ci     <- c(50,100,150,200,300,500,800,1000,1200,1500)
  PPFDs  <- c(50, 150, 300, 500, 1200,1800)

  # Convert into a list of vectors (like a data frame but plain list)

  #datatype == "ACi(different Q)"
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
  return(list(Ci=df$C_i, PPFD=df$PPFD, An=An, Phi2=Phi2))
}


add_noise <- function(An,sd_abs=0.2, sd_rel=0.05){

  A_sd <- sqrt(sd_abs^2 + (sd_rel * An)^2)
  A_obs <- pmax(An + rnorm(length(An), 0, A_sd), -2)
  return(A_obs)
}


datatypes <- c("ACi(satQ)", "ACi & AQ", "ACi(different Q)", "ACi(different Q) & CF")


row_id  <- 1
simulated_curves <- list()
for (rep_x in 1:100){
  print(paste0("Iteration: ", rep_x))
  true_pars <- sampled_params[[rep_x]]
  for (m in 1:length(rst_functions)) {
    model_m <- rst_functions[[m]]

    sim_res0 <- generate_curves(model_m,true_pars,names(rst_functions)[m])
    noisy_An <- add_noise(sim_res0$An)
    noisy_Phi2 <- add_noise(sim_res0$Phi2,sd_abs = 0.05, sd_rel = 0.01)
    if (m==9){
      dl <- 4
    }else{dl <- 3}
    for (d in 1:dl){

      if (d==1){ # datatype=="ACi(satQ)"
        sim_res <- list(Ci=sim_res0$Ci[sim_res0$PPFD==1800],
                        PPFD=sim_res0$PPFD[sim_res0$PPFD==1800],
                        An=noisy_An[sim_res0$PPFD==1800])

      }else if(d==2){# datatype == "ACi & AQ"
        sim_res <- list(Ci=c(sim_res0$Ci[sim_res0$PPFD==1800],sim_res0$Ci[sim_res0$Ci==300]),
                        PPFD=c(sim_res0$PPFD[sim_res0$PPFD==1800],sim_res0$PPFD[sim_res0$Ci==300]),
                        An=c(noisy_An[sim_res0$PPFD==1800],noisy_An[sim_res0$Ci==300]))

      }else if(d==3){ #"ACi(different Q)"
        sim_res <- list(Ci=sim_res0$Ci,
                        PPFD=sim_res0$PPFD,
                        An=noisy_An)
      }else if (d==4){ #"ACi(different Q) & CF"
        sim_res <- sim_res0
        sim_res$An <- noisy_An
        sim_res$Phi2 <- noisy_Phi2
      }

      simulated_curves[[row_id]] <- tibble(
        Ci       = sim_res$Ci,
        PPFD     = sim_res$PPFD,
        An       = sim_res$An,
        Phi2     = sim_res$Phi2
      )

      names(simulated_curves)[row_id] <- paste0("rep",rep_x,"_",names(rst_functions)[m],"_",datatypes[d])
      row_id <- row_id + 1

    }
  }
}

save(simulated_curves,file="results/Rdata/fig2_simulated_curves.RData")
