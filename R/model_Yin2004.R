#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' Wp: TPU-limited carboxylation rate
#' Reference: Extension of a biochemical model for the generalized  stoichiometry of electron transport limited C3  photosynthesis
#' @param envs List of environmental variables:
#'   \describe{
#'     \item{`PPFD`}{Photosynthetic photon flux density in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`C_i`}{CO2 concentration in intra-cellular space in μmol mol\eqn{^{-1}}.}
#'   }
#' @param pars List of parameters:
#'   \describe{
#'     \item{`V_cmax`}{Maximum rate of RuBP-saturated carboxylation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`J_max`}{Maximum rate of electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as O2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
#'     \item{`f_cyc`}{Fraction of cyclic electron transport around PSI, unitless.}
#'     \item{`f_Q`}{Fraction of electron transport that follows the Q-cycle, unitless.})
#'     \item{`f_pseudo`}{Fraction of pseudocyclic electron transport, unitless.})
#'     \item{`phi2m`}{Maximum electron transport efficiency of PSII in mol e– mol\eqn{^{-1}} PPFD.}
#'   }
#' @return List of calculated rate:
#'   \describe{
#'     \item{`W_c`}{RuBP-saturated carboxylation rate in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`W_j`}{RuBP regeneration-limited carboxylation rate in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`W_p`}{TPU-limited carboxylation rate in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Wmin`}{Minimum between W_c, W_j, and W_p.}
#'     \item{`An`}{Net assimilation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Ac`}{Net assimilation, RuBP-saturated, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Aj`}{Net assimilation, RuBP regeneration-limited, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Ap`}{Net assimilation, TPU-limited, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Limitation`}{Character vector specifying the limitation at which CO2 level.}
#'   }
#' @export


Yin2004 = function(envs,pars) {

  a <- 4 * envs$C_i + 8 * pars$gamma_star
  b <- 3 * envs$C_i + 7 * pars$gamma_star

  # calculate one of the fractions, with two others given
  # eq 1)
  if (is.na(pars$f_cyc)){
    pars$f_cyc <- (a*(2 + pars$f_Q)/(pars$h*b) - (1 - pars$f_pseudo)) / (a/(pars$h*b) - 1)
  }else if (is.na(pars$f_Q)){
    pars$f_Q <- pars$h*b*(1 - pars$f_cyc - pars$f_pseudo)/a - 2 + pars$f_cyc
  }else if (is.na(pars$f_pseudo)){
    pars$f_pseudo <- 1 - pars$f_cyc - a*(2 + pars$f_pseudo - pars$f_cyc)/(h*b)
  }

  alpha_J <- (1 - pars$f_cyc) / (1 + ((1 - pars$f_cyc) / pars$phi2m))

  B <- - (alpha_J * envs$PPFD + pars$J_max)
  J <- (- B - sqrt((B)^2 - 4 * pars$theta_J * alpha_J * envs$PPFD * pars$J_max)) / (2 * pars$theta_J)

  pars$K_CO <- pars$K_C * (envs$O / pars$K_O)
  ret = list(
    Wc = pars$V_cmax * envs$C_i / (envs$C_i + pars$K_CO),
    Wj = J * ((2 + pars$f_Q - pars$f_cyc) * envs$C_i) / (pars$h * (3 * envs$C_i + 7 * pars$gamma_star) * (1 - pars$f_cyc))
  )

  Wmin <- c()
  limiting_factor <- character(length(envs$C_i))
  for (i in seq_along(envs$C_i)) {
    Wmin[i] <- min(ret$Wc[i], ret$Wj[i])
    limiting_factor[i] <- c("Wc", "Wj")[which.min(c(ret$Wc[i], ret$Wj[i]))]
  }

  ret$Wmin = Wmin
  ret$An = (1- pars$gamma_star/envs$C_i)*Wmin - pars$R_d

  ret$Ac = (1- pars$gamma_star/envs$C_i)*ret$Wc - pars$R_d
  ret$Aj = (1- pars$gamma_star/envs$C_i)*ret$Wj - pars$R_d

  ret$Limitation <- limiting_factor

  return(ret)
}
