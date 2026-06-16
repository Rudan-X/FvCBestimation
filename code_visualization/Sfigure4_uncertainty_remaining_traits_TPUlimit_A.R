library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)


devtools::load_all()


load(file="results/Rdata/fig2_sampled_params.RData")


call_function_FvCB80 <- function(Ci, O, PPFD, pars) {
  model_FvCB80(Ci, 0L, 0, O, PPFD,
               pars$V_cmax, pars$J_max, pars$R_d,
               pars$K_C, pars$K_O, pars$gamma_star,
               pars$alpha_J, pars$theta_J)
}

call_function_Harley92 <- function(Ci, O, PPFD, pars) {
  model_Harley92(Ci, O, PPFD,
                 pars$V_cmax, pars$J_max, pars$R_d,
                 pars$K_C, pars$K_O, pars$gamma_star,
                 pars$alpha_J,pars$gm)
}

call_function_Caemmerer00 <- function(Ci, O, PPFD, pars) {
  model_Caemmerer00(Ci, 0L, 0, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
                    pars$alpha_J, pars$theta_J,
                    pars$V_tpu, pars$alpha_G)
}

call_function_Ethier04 <- function(Ci, O, PPFD, pars) {
  model_Ethier04(Ci, O, PPFD,
                 pars$V_cmax, pars$J_max, pars$R_d,
                 pars$K_C, pars$K_O, pars$gamma_star,
                 pars$alpha_J, pars$theta_J,
                 pars$gm)
}

call_function_Yin04 <- function(Ci, O, PPFD, pars) {
  model_Yin04(Ci, 0L, 0, O, PPFD,
              pars$V_cmax, pars$J_max, pars$R_d,
              pars$K_C, pars$K_O, pars$gamma_star,
              pars$phi2m, pars$theta_J,
              pars$f_Q, pars$f_pseudo, pars$h)
}


call_function_Dubois07 <- function(Ci, O, PPFD, pars) {
  model_Dubois07(Ci, 0L, 0, O, PPFD,
                 pars$V_cmax, pars$J_max, pars$R_d,
                 pars$K_C, pars$K_O, pars$gamma_star,
                 pars$alpha_J, pars$theta_J,
                 pars$V_tpu, pars$alpha_G)
}


call_function_Tholen12 <- function(Ci, O, PPFD, pars) {
  model_Tholen12(Ci, O, PPFD,
                 pars$V_cmax, pars$J_max, pars$R_d,
                 pars$K_C, pars$K_O, pars$gamma_star,
                 pars$alpha_J,pars$theta_J, pars$g_ch ,pars$g_wp)
}


call_function_Busch18 <- function(Ci, O, PPFD, pars) {
  model_Busch18(Ci, 0L, 0, O, PPFD,
                pars$V_cmax, pars$J_max, pars$R_d,
                pars$K_C, pars$K_O, pars$gamma_star,
                pars$alpha_J, pars$theta_J,
                pars$V_tpu, pars$N_max, pars$max_aG, pars$max_aS)
}

call_function_Xiao20 <- function(Ci, O, PPFD, pars) {
  model_Xiao20(Ci, O, PPFD,
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

  out <- t(vapply(seq_along(envs$C_i), function(i) {
    r <- model(envs$C_i[i], envs$O[i], envs$PPFD[i], parameters)
    c(An = as.character(r$An), Limitation = r$Limitation)
  }, character(2)))

  out <- as.data.frame(out)
  out$An <- as.numeric(out$An)


  if (model_name == "Xiao20"){

    Phi2 <- vapply(seq_along(envs$C_i), function(i)
      model(envs$C_i[i], envs$O[i], envs$PPFD[i],
            parameters)$Phi2,
      numeric(1)
    )
  }else{Phi2=NA_real_}
  return(list(Ci=df$C_i, PPFD=df$PPFD, An=out$An, Phi2=Phi2, Limitation=out$Limitation))
}


# add_noise <- function(An,sd_abs=0.2, sd_rel=0.05){
#
#   A_sd <- sqrt(sd_abs^2 + (sd_rel * An)^2)
#   A_obs <- pmax(An + rnorm(length(An), 0, A_sd), -2)
#   return(A_obs)
# }


datatypes <- c("ACi(satQ)", "ACi & AQ", "ACi(different Q)", "ACi(different Q) & CF")


row_id  <- 1
simulated_curves <- list()
for (rep_x in 1:100){
  print(paste0("Iteration: ", rep_x))
  true_pars <- sampled_params[[rep_x]]
  for (m in 1:length(rst_functions)) {
    model_m <- rst_functions[[m]]

    sim_res0 <- generate_curves(model_m,true_pars,names(rst_functions)[m])
    # noisy_An <- add_noise(sim_res0$An)
    # noisy_Phi2 <- add_noise(sim_res0$Phi2,sd_abs = 0.05, sd_rel = 0.01)
    if (m==9){
      dl <- 4
    }else{dl <- 3}
    for (d in 1:dl){

      if (d==1){ # datatype=="ACi(satQ)"
        sim_res <- list(Ci=sim_res0$Ci[sim_res0$PPFD==1800],
                        PPFD=sim_res0$PPFD[sim_res0$PPFD==1800],
                        An=sim_res0$An[sim_res0$PPFD==1800],
                        Limitation=sim_res0$Limitation[sim_res0$PPFD==1800])

      }else if(d==2){# datatype == "ACi & AQ"
        sim_res <- list(Ci=c(sim_res0$Ci[sim_res0$PPFD==1800],sim_res0$Ci[sim_res0$Ci==300]),
                        PPFD=c(sim_res0$PPFD[sim_res0$PPFD==1800],sim_res0$PPFD[sim_res0$Ci==300]),
                        An=c(sim_res0$An[sim_res0$PPFD==1800],sim_res0$An[sim_res0$Ci==300]),
                        Limitation=c(sim_res0$Limitation[sim_res0$PPFD==1800],sim_res0$Limitation[sim_res0$Ci==300]))

      }else if(d==3){ #"ACi(different Q)"
        sim_res <- list(Ci=sim_res0$Ci,
                        PPFD=sim_res0$PPFD,
                        An=sim_res0$An,
                        Limitation=sim_res0$Limitation)
      }else if (d==4){ #"ACi(different Q) & CF"
        sim_res <- sim_res0
        sim_res$An <- sim_res0$An
        sim_res$Phi2 <- sim_res0$Phi2
      }

      simulated_curves[[row_id]] <- tibble(
        Ci       = sim_res$Ci,
        PPFD     = sim_res$PPFD,
        An       = sim_res$An,
        Phi2     = sim_res$Phi2,
        Limitation = sim_res$Limitation
      )

      names(simulated_curves)[row_id] <- paste0("rep",rep_x,"_",names(rst_functions)[m],"_",datatypes[d])
      row_id <- row_id + 1

    }
  }
}


check <- rep(NA,length(simulated_curves))
for (i in 1:length(check)){
  lim <- simulated_curves[[i]]$Limitation
  check[i] <- sum(grepl("p", lim))
}

ind <- which(check!=0)
p_sim <- simulated_curves[ind]
p_sim <- p_sim[grepl("ACi\\(satQ\\)", names(p_sim))]

tpu_ind <- unique(sub(".*(rep\\d+).*", "\\1", names(p_sim)))
tpu_ind <- as.numeric(gsub("rep","",tpu_ind))

save(tpu_ind, file="results/Rdata/data_visualization/Sfig_TPU_ind.RData")

library(dplyr)
library(purrr)
library(ggplot2)

df <- imap_dfr(p_sim, ~ mutate(.x, key = .y))
df$model <- sub("^rep\\d+_([^_]+)_.*$", "\\1", df$key)

df$model[df$model== "Caemmerer00"] <-  "vonC00"
# df$model[df$model== "Busch18"] <-  "Busch17"
df$model <- factor(df$model,levels=c("vonC00", "Dubois07","Busch18"))
df$Limitation <- factor(df$Limitation,
  levels = c( "Wc", "Wj", "Wp","Ac", "Aj", "Ap")
)
# g1 <- ggplot(df, aes(x = Ci, y = An,  group = key, shape = Limitation)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.5) +
#   theme_bw() +
#   facet_grid(model~.)+
#   labs(x = "Ci", y = "An") +
#   theme(legend.position = "bottom")


g1 <- ggplot(df, aes(x = Ci, y = An,
               color = key,
               shape = Limitation)) +
  geom_line(aes(group = key), linewidth = 0.5) +
  geom_point(size = 2) +
  facet_grid(.~model)+
  theme_bw()+
  scale_shape_manual(
    values = c(Wc = 16, Wj = 17, Wp = 15, Ac = 1, Aj = 2, Ap = 0)
  ) +
  theme(legend.position = "bottom") +
  labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
       y = expression(An~"(µmol m"^{-2}~s^{-1}*")"),
       color = "", shape="") +
  guides(color = "none", shape = guide_legend(nrow = 1, byrow = TRUE))

g1
