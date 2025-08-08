#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' Reference: Theoretical Considerations when Estimating the Mesophyll Conductance to CO2 Flux by Analysis of the Response of Photosynthesis to CO2
#' @param envs List of environmental variables:
#'   \describe{
#'     \item{`PPFD`}{Photosynthetic photon flux density in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`O`}{O2 concentration in chloroplast in μmol mol\eqn{^{-1}}.}
#'     \item{`C_i`}{CO2 concentration in intracellular space in μmol mol\eqn{^{-1}}.}
#'
#'   }
#' @param pars List of parameters:
#'   \describe{
#'     \item{`V_cmax`}{Maximum rate of RuBP-saturated carboxylation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`J_max`}{Maximum rate of electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as CO2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`gm`}{Mesophyll conductance in mol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'   }
#' @return List of calculated rate:
#'   \describe{
#'     \item{`An`}{Net assimilation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Ac`}{Net assimilation, RuBP-saturated, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Aj`}{Net assimilation, RuBP regeneration-limited, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Limitation`}{Character vector specifying the limitation at which CO2 level.}
#'   }
#' @export


Harley1992 = function(envs,pars) {
  ret = list()
  Amin <- numeric(length(envs$C_i))
  limiting_factor <- character(length(envs$C_i))
  ret$Ac <- numeric(length(envs$C_i))
  ret$Aj <- numeric(length(envs$C_i))

  J = 0.24 * envs$PPFD / sqrt(1 + (0.24 * envs$PPFD/pars$J_max)^2)

  for (i in 1:length(envs$C_i)) {
    Cc_guess <- envs$C_i[i]
    continue <- TRUE
    iter <- 0
    while (continue & iter<1000){
      Ac <- pars$V_cmax * (Cc_guess - pars$gamma_star) / (Cc_guess + (pars$K_C * (1 + envs$O[i] / pars$K_O)))
      Aj <- J[i] * (Cc_guess - pars$gamma_star) / (4 * (Cc_guess + 2 * pars$gamma_star))
      A_new <- min(Ac, Aj) - pars$R_d
      Cc_new <- envs$C_i[i] - A_new / pars$gm
      if (abs(Cc_new - Cc_guess) < 0.001){
        continue <- FALSE
      }else{
        Cc_guess <- Cc_new
        iter <- iter+1
      }
    }

    if (iter==999){
      Amin[i] <- NA_real_
      ret$Ac[i] <- NA_real_
      ret$Aj[i] <- NA_real_
      limiting_factor[i] <- "Failed"

    }else{
      Amin[i] <- A_new
      ret$Ac[i] <- Ac
      ret$Aj[i] <- Aj
      limiting_factor[i] <- c("Ac", "Aj")[which.min(c(ret$Ac[i], ret$Aj[i]))]
    }

  }

  ret$An = Amin
  ret$Limitation <- limiting_factor

  return(ret)
}

