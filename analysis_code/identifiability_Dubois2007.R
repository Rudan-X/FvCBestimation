# https://rpubs.com/DeTwes/Identifiability-Analysis
library(tidyr)
library(dplyr)
library(ggpubr)

devtools::load_all()


calc_Anet <- function(Ci, Vcmax, J, Rd, gamma_s){ #, KC, KO
  envs <- list(
    C_chl = Ci,
    O = rep(210, length(Ci)))

  pars <- list(
    R_d= Rd,
    K_C= 400,
    K_O= 280,
    gamma_star= gamma_s,
    V_cmax= Vcmax,
    J = J
  )

  C_transition=NULL
  Ares <- Dubois2007(envs=envs,pars=pars,Ct=C_transition)
  return(Ares$An)
}


C <- c(50,100,150,200,300,400,500,600,700,800,900,1000, 1200)
true.val <- c(150, 200, 1, 45) #, 400, 280
An <- calc_Anet(C, true.val[1], true.val[2], true.val[3], true.val[4]) #, true.val[4], true.val[5]

ACi_curve <- data.frame(Ci = C, A = An)

plot(C, An)

start_vals <- list(Vcmax = 50, J = 100, Rd = 0.1, gamma_s = 20) # , KO = 400, KC = 280

lb0 <- c(Vcmax = 10, J = 10, Rd = 0, gamma_s = 10 ) #, KO = 10, KO = 10
ub0 <- c(Vcmax = 500, J = 500, Rd = 50, gamma_s = 100) #, KO = 700, KO = 700
fit <- nls(
  A ~ calc_Anet(Ci, Vcmax, J, Rd, gamma_s), #, KO, KC
  data = ACi_curve,
  start = start_vals,
  algorithm = "port",
  lower = lb0,
  upper = ub0,
  control = nls.control(maxiter = 500)
)

best <- coef(fit)


args <- list(
  Ci     = ACi_curve$Ci,
  Vcmax  = best["Vcmax"],
  J      = best["J"],
  Rd     = best["Rd"],
  gamma  = best["gamma_s"]
)

# Parameters to profile
params <- c("Vcmax", "J", "Rd", "gamma_s")
profiles <- list()


for (p in params) {

  # Create a sequence of fixed values for the parameter
  p_vals <- seq(best[p] * 0.5, best[p] * 1.5, length.out = 20)

  # Store log-likelihood values
  logL_vals <- numeric(length(p_vals))
  re_opt <- setdiff( params, p)
  start <- best[re_opt]
  lb <- lb0[re_opt]
  ub <- ub0[re_opt]
  for (i in seq_along(p_vals)) {
    try({

      param_exprs <- lapply(params, function(param) {
        if (param == p) p_vals[i] else as.name(param)
      })
      names(param_exprs) <- params

      formula <- as.formula(
        substitute(
          A ~ calc_Anet(Ci, Vcmax, J, Rd, gamma_s),
          param_exprs
        )
      )

      fit_tmp <- nls(formula = formula,
        data = ACi_curve,
        start = start,
        algorithm = "port",
        lower = lb,
        upper = ub,
        control = nls.control(maxiter = 500)
      )

      rss <- sum(resid(fit_tmp)^2)
      sigma2 <- rss / nrow(ACi_curve)
      logL_vals[i] <- -nrow(ACi_curve) / 2 * log(sigma2)
    })
  }

  # Save results for plotting
  profiles[[p]] <- data.frame(
    param = p,
    fixed_value = p_vals,
    logLik = logL_vals
  )
}


profiles_df <- do.call(rbind, profiles)

# Plot profiles
library(ggplot2)
ggplot(profiles_df, aes(x = fixed_value, y = logLik)) +
  geom_line() +
  geom_point() +
  facet_wrap(~param, scales = "free_x") +
  geom_vline(data = data.frame(param = params, best = best),
             aes(xintercept = best), col = "red", linetype = 2) +
  theme_minimal() +
  labs(x = "Fixed parameter value", y = "Log-likelihood",
       title = "Identifiability Profiles")
