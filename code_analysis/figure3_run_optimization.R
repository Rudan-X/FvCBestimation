
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

devtools::load_all()

FvCB80 <- model_FvCB1980()
Caemmerer00 <- model_Caemmerer2000(base = FvCB80, ci_tpu_min = 100)
Ethier04 <- model_Ethier2004()

models <- list(FvCB80 = FvCB80,
               Harley92 = model_Harley1992(base = FvCB80),
               Caemmerer00  = Caemmerer00,
               Ethier04  = Ethier04,
               Yin04 = model_Yin2004(base = FvCB80),
               Dubois07 = model_Dubois2007(base = Caemmerer00, ci_tpu_min = 100),
               Tholen12 = model_Tholen2012(base = FvCB80),
               Busch18 = model_Busch2018(base = FvCB80, ci_tpu_min = 100),
               Xiao20 = model_Xiao2021(base = Ethier04))

fit_names <-c("V_cmax","J_max","R_d","gamma_star","theta_J","alpha_J", "K_C", "K_O")

model_params <- list(FvCB80 = fit_names,
               Harley92 = c(fit_names, "gm"),
               Caemmerer00  = c(fit_names, "V_tpu", "alpha_G"),
               Ethier04  = c(fit_names, "gm"),
               Yin04 = c(fit_names, "h", "phi2m","f_pseudo","f_Q"),
               Dubois07 = c(fit_names, "V_tpu", "alpha_G"),
               Tholen12 = c(fit_names, "g_ch", "g_wp"),
               Busch18 = c(fit_names, "max_aG", "max_aS","N_max"),
               Xiao20 = c(fit_names, "gm", "s","Phi2LL"))

# 1) Environmental condition to simulate A

Ci     <- c(50,100,150,200,300,500,600,700,800,1000,1200)
PPFDs  <- c(50, 100, 500,1500)
envs_template <- list(C_i = Ci, O = rep(210, length(Ci)))


# 2) Reasonable starts / bounds FOR ALL POSSIBLE PARAMETERS
#    (models will pick the subset they need)

lower_all <- c(
  V_cmax=50, J_max=100, R_d=0.5, gamma_star=10, V_tpu=8,
  gm=0.1, alpha_J=0.3,  theta_J=0.5, K_C=200, K_O=150, alpha_G=0.3,
  phi2m = 0.5,  f_Q =  0.5,  f_pseudo = 0,
  h = 3,  max_aG  = 0.01,  max_aS =0.01,  N_max = 0.1,
  g_ch = 0.05,  g_wp =0.1,  s = 0.6)



upper_all <- c(
  V_cmax=200, J_max=250, R_d=5, gamma_star=60,V_tpu=20,
  gm=0.7, alpha_J=0.5,  theta_J=0.9, K_C=300, K_O=300, alpha_G=0.8,
  phi2m = 0.9, f_Q =  1, f_pseudo = 0.3,
  h = 14/3, max_aG  =  0.2, max_aS =0.5, N_max = 2, g_ch = 0.7,
  g_wp = 2, s = 0.6)


sample_true_pars <- function(){
  temp <- list(
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

  temp$Phi2LL <- temp$alpha_J / temp$s
  return(temp)
}

load(file="results/Rdata/fig2_sampled_params.RData")
load(file="results/Rdata/fig2_simulated_curves.RData")

for (rep_x in 1:100){

  for (i in 1:length(res_list)){
    case <- names(simulated_curves)[i]
    check <- strsplit(case,"_")
    datatype <- check[[1]][3]
    model <- check[[1]][2]

    true_curve <- simulated_curves[[case]]

    envs <- list(C_i= true_curve$Ci, PPFD= true_curve$PPFD, O= rep(210, nrow(true_curve)))

  }
}

fit_one_repitition <- function(rep_id,models,   # list of compiled Stan models
                               iter = 3000, warmup = 1000, chains = 4,
                               adapt_delta = 0.95, max_treedepth = 12) {

  res_list <- list()
  # storage
  for (m in 1:length(models)) {
    model_m <- models[[m]]
    sampled_start <- sample_true_pars() # sample random starting point
    start <- unlist(sampled_start[fit_names])
    lower <- lower_all[fit_names]
    upper <- upper_all[fit_names]

    full_par <- sample_true_pars()
    full_par <- full_par[model_params[[m]]]

    true_curve <- simulated_curves[[case]]

    rss_fn <- function(theta_vec){
      full_par[names(theta_vec)] <- as.list(theta_vec)
      preds <- predict_A(model_m, envs, full_par)$An

      return(sum((true_curve$An - preds)^2))
    }
    opt <- optim(par = start, fn = rss_fn, method = "L-BFGS-B",
                 lower = lower, upper = upper, control = list(maxit = 2000))

    opt$value
    # helper: vectorized prediction per PPFD (no per-row predict_A)
    predict_all <- function(pars_try){
      preds <- numeric(nrow(dat))
      for (pp in unique(dat$PPFD)) {
        idx  <- which(dat$PPFD == pp)
        envs <- list(C_i = dat$Ci[idx], O = rep(210, length(idx)), PPFD = rep(pp, length(idx)))
        val  <- try(predict_A(mdl_fun, envs, pars_try)$An, silent = TRUE)
        if (inherits(val, "try-error") || any(!is.finite(val)) || length(val) != length(idx))
          return(rep(NaN, nrow(dat)))
        preds[idx] <- as.numeric(val)
      }
      return(preds)
    }

    # objective
    rss_local <- function(theta){
      input_pars <- modifyList(full_pars, as.list(theta))

      preds <- predict_all(input_pars)
      if (any(!is.finite(preds))){
        return(1e12)
      }
      rss <- sum((dat$A_obs  - preds)^2)
      if (!is.finite(rss)){
        return(1e12)
      }else{
        return(rss)
      }
    }
  }
}
