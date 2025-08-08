#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' On the need to incorporate sensitivity to CO2 transfer  conductance into the Farquhar–von Caemmerer–Berry leaf  photosynthesis model
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
#'     \item{`gamma_star`}{Chloroplastic CO2 photocompensation point, same units as CO2.}
#'     \item{`alpha_J`}{quantum efficiency, unitless.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
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


EthierLong2004 = function(envs,pars) {
  ret = list()
  Amin <- numeric(length(envs$C_i))
  limiting_factor <- character(length(envs$C_i))
  ret$Ac <- numeric(length(envs$C_i))
  ret$Aj <- numeric(length(envs$C_i))

  pars$K_CO <- pars$K_C * (1 + envs$O / pars$K_O)

  b <- - (pars$alpha_J * envs$PPFD + pars$J_max)
  J <- (- b - sqrt((b)^2 - 4 * pars$theta_J * pars$alpha_J * envs$PPFD * pars$J_max)) / (2 * pars$theta_J)


  # quadratic form of Ac from equation 3
  a <- -1 / pars$gm
  b <- ((pars$V_cmax - pars$R_d) / pars$gm) + envs$C_i + pars$K_CO
  c <- pars$R_d * (envs$C_i + pars$K_CO) - pars$V_cmax * (envs$C_i - pars$gamma_star)

  sqrt_term <- sqrt(b^2 - 4 * a * c)
  ret$Ac <- (-b + sqrt_term) / (2 * a)

  # quadratic form of Aj from equation
  a <- -4 / pars$gm
  b <- 4 * (envs$C_i + 2 * pars$gamma_star) - 4 * pars$R_d/pars$gm +  J/pars$gm
  c <- 4 * pars$R_d * (envs$C_i + 2 * pars$gamma_star) - J * (envs$C_i - pars$gamma_star)

  sqrt_term <- sqrt(b^2 - 4 * a * c)
  ret$Aj <- (-b + sqrt_term) / (2 * a)


  for (i in seq_along(envs$C_i)) {
    Amin[i] <- min(ret$Ac[i], ret$Aj[i])
    limiting_factor[i] <- c("Ac", "Aj")[which.min(c(ret$Ac[i], ret$Aj[i]))]
  }

  ret$An = Amin
  ret$Limitation <- limiting_factor

  return(ret)
}

