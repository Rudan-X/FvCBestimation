#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' Reference: Optimizing the statistical estimation of the parameters of  the Farquhar–von Caemmerer–Berry model of  photosynthesis
#' @param envs List of environmental variables:
#'   \describe{
#'     \item{`PPFD`}{Photosynthetic photon flux density in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`C_i`}{Intercellular CO2 concentration in μmol mol\eqn{^{-1}}.}
#'     \item{`O`}{O2 concentration in chloroplast in μmol mol\eqn{^{-1}}.}
#'   }
#' @param pars List of parameters:
#'   \describe{
#'     \item{`V_cmax`}{Maximum rate of RuBP-saturated carboxylation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`J_max`}{Maximum rate of electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as O2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`alpha_J`}{quantum efficiency, unitless.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
#'     \item{`V_tpu`}{Rate of triose phosphate export from the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`alpha_tpu`}{Fraction of  the glycolate carbon that is not returned to the chloroplast, unitless.}
#'
#'   }
#' @param Ct Optional input (default = `NULL`). A numeric vector of two element:
#' first is Ci transition between the Rubisco and RuBP Regeneration limitations
#' second is Ci transition between the RuBP Regeneration and TPU limitations
#' @return List of calculated rate:
#'   \describe{
#'     \item{`An`}{Net assimilation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Ac`}{Net assimilation, RuBP-saturated, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Aj`}{Net assimilation, RuBP regeneration-limited, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Limitation`}{Character vector specifying the limitation at which CO2 level.}
#'   }
#' @export


Dubois2007 = function(envs,pars,Ct=NULL) {

  bJ <- - (pars$alpha_J * envs$PPFD + pars$J_max)
  J <- (- bJ - sqrt((bJ)^2 - 4 * pars$theta_J * pars$alpha_J * envs$PPFD * pars$J_max)) / (2 * pars$theta_J)

  pars$K_CO <- pars$K_C * (1 + envs$O / pars$K_O)
  ret = list(
    Wc = pars$V_cmax * envs$C_i / (envs$C_i + pars$K_CO),
    Wj = J * envs$C_i / (4 * envs$C_i + 8 * pars$gamma_star),
    # Wp = 3 * pars$V_tpu * envs$C_i / (envs$C_i - pars$gamma_star)
    Wp = 3 * pars$V_tpu * envs$C_i / (envs$C_i - (1 + 3*pars$alpha_tpu) * pars$gamma_star)
  )

  ret$Ac = (1- pars$gamma_star/envs$C_i)*ret$Wc - pars$R_d
  ret$Aj = (1- pars$gamma_star/envs$C_i)*ret$Wj - pars$R_d
  ret$Ap = (1- pars$gamma_star/envs$C_i)*ret$Wp - pars$R_d

  ret$An <- numeric(length(envs$C_i))
  ret$Limitation <- character(length(envs$C_i))

  if (is.null(Ct)){
    for (i in seq_along(envs$C_i)) {
      if (envs$C_i[i] > pars$gamma_star) {
        ret$An[i] <- min(ret$Ac[i], ret$Aj[i],ret$Ap[i])
        ret$Limitation[i] <- c("Ac", "Aj","Ap")[which.min(c(ret$Ac[i], ret$Aj[i], ret$Ap[i]))]
      } else {
        ret$An[i] <- min(ret$Ac[i], ret$Aj[i])
        ret$Limitation[i] <- c("Ac", "Aj")[which.min(c(ret$Ac[i], ret$Aj[i]))]
      }

    }
  }else{
    Ci_cj <- Ct[1]# 450
    Ci_jp <- Ct[2]# 700

    for (i in seq_along(envs$C_i)) {
      if (envs$C_i[i] < Ci_cj ) {
        ret$An[i] <- ret$Ac[i]
        ret$Limitation[i] <- "Ac"
      } else if (envs$C_i[i] >= Ci_cj & envs$C_i[i] <= Ci_jp ) {
        ret$An[i] <- ret$Aj[i]
        ret$Limitation[i] <- "Aj"
      } else {
        ret$An[i] <- ret$Aj[i]
        ret$Limitation[i] <- "Ap"
      }
    }
  }

  return(ret)
}

