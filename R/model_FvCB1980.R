#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' Wp: TPU-limited carboxylation rate
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
#'     \item{`V_tpu`}{Rate of triose phosphate export from the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as O2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`alpha_J`}{quantum efficiency, unitless.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
#'   }
#' @param Ct Optional input (default = `NULL`). A numeric vector of two element:
#' first is Ci transition between the Rubisco and RuBP Regeneration limitations
#' second is Ci transition between the RuBP Regeneration and TPU limitations
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


FvCB1980 = function(envs,pars,Ct=NULL) {

  bJ <- - (pars$alpha_J * envs$PPFD + pars$J_max)
  J <- (- bJ - sqrt((bJ)^2 - 4 * pars$theta_J * pars$alpha_J * envs$PPFD * pars$J_max)) / (2 * pars$theta_J)


  K_CO <- pars$K_C * (envs$O / pars$K_O)
  ret = list(
    Wc = pars$V_cmax * envs$C_i / (envs$C_i + K_CO),
    Wj = J * envs$C_i / (4 * envs$C_i + 8 * pars$gamma_star),
    Wp = 3 * pars$V_tpu * envs$C_i / (envs$C_i - pars$gamma_star)
  )

  Wmin <- c()
  limiting_factor <- character(length(envs$C_i))

  if (is.null(Ct)){
    for (i in seq_along(envs$C_i)) {
      if (envs$C_i[i] > pars$gamma_star) {
        Wmin[i] <- min(ret$Wc[i], ret$Wj[i], ret$Wp[i])
        limiting_factor[i] <- c("Wc", "Wj", "Wp")[which.min(c(ret$Wc[i], ret$Wj[i], ret$Wp[i]))]
      } else {
        Wmin[i] <- min(ret$Wc[i], ret$Wj[i])
        limiting_factor[i] <- c("Wc", "Wj")[which.min(c(ret$Wc[i], ret$Wj[i]))]
      }
    }
  }else{
    Ci_cj <- Ct[1]# 450
    Ci_jp <- Ct[2]# 700

    for (i in seq_along(envs$C_i)) {
      if (envs$C_i[i] < Ci_cj ) {
        Wmin[i] <- ret$Wc[i]
        limiting_factor[i] <- "Wc"
      } else if (envs$C_i[i] >= Ci_cj & envs$C_i[i] <= Ci_jp ) {
        Wmin[i] <- ret$Wj[i]
        limiting_factor[i] <- "Wj"
      }else{
        Wmin[i] <- ret$Wp[i]
        limiting_factor[i] <- "Wp"
      }
    }
  }

  ret$Wmin = Wmin
  ret$An = (1- pars$gamma_star/envs$C_i)*Wmin - pars$R_d

  ret$Ac = (1- pars$gamma_star/envs$C_i)*ret$Wc - pars$R_d
  ret$Aj = (1- pars$gamma_star/envs$C_i)*ret$Wj - pars$R_d
  ret$Ap = (1- pars$gamma_star/envs$C_i)*ret$Wp - pars$R_d
  ret$Limitation <- limiting_factor

  return(ret)
}


