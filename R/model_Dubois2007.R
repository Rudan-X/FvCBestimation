#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' @param envs List of environmental variables:
#'   \describe{
#'     \item{`C_chl`}{CO2 concentration in chloroplast in μmol mol\eqn{^{-1}}.}
#'     \item{`O`}{O2 concentration in chloroplast in μmol mol\eqn{^{-1}}.}
#'   }
#' @param pars List of parameters:
#'   \describe{
#'     \item{`V_cmax`}{Maximum rate of RuBP-saturated carboxylation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`J`}{Electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as CO2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'   }
#' @param Ct Optional input (default = `NULL`). A numeric vector of two element:
#' first is Ci transition between the Rubisco and RuBP Regeneration limitations
#' second is Ci transition between the RuBP Regeneration and TPU limitations
#' @return List of calculated rate:
#'   \describe{
#'     \item{`W_c`}{RuBP-saturated carboxylation rate in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`W_j`}{RuBP regeneration-limited carboxylation rate in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Wmin`}{Minimum between W_c, W_j, and W_p.}
#'     \item{`An`}{Net assimilation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Ac`}{Net assimilation, RuBP-saturated, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Aj`}{Net assimilation, RuBP regeneration-limited, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Limitation`}{Character vector specifying the limitation at which CO2 level.}
#'   }
#' @export


Dubois2007 = function(envs,pars,Ct=NULL) {

  ret = list(
    Wc = pars$V_cmax * envs$C_chl / (envs$C_chl + pars$K_C * (1 + envs$O / pars$K_O)),
    Wj = pars$J * envs$C_chl / (4 * envs$C_chl + 8 * pars$gamma_star)
  )

  ret$Ac = (1- pars$gamma_star/envs$C_chl)*ret$Wc - pars$R_d
  ret$Aj = (1- pars$gamma_star/envs$C_chl)*ret$Wj - pars$R_d

  Amin <- c()
  limiting_factor <- character(length(envs$C_chl))

  if (is.null(Ct)){
    for (i in seq_along(envs$C_chl)) {
      Amin[i] <- min(ret$Ac[i], ret$Aj[i])
      limiting_factor[i] <- c("Ac", "Aj")[which.min(c(ret$Ac[i], ret$Aj[i]))]
    }
  }else{
    Ci_cj <- Ct[1]# 450
    Ci_jp <- Ct[2]# 700

    for (i in seq_along(envs$C_chl)) {
      if (envs$C_chl[i] < Ci_cj ) {
        Amin[i] <- ret$Ac[i]
        limiting_factor[i] <- "Ac"
      } else {
        Wmin[i] <- ret$Aj[i]
        limiting_factor[i] <- "Aj"
      }
    }
  }

  # ret$Wmin = Wmin
  ret$An = Amin
  ret$Limitation <- limiting_factor

  return(ret)
}

