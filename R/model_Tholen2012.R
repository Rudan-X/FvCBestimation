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
#'     \item{`S_co`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`alpha_J`}{quantum efficiency, unitless.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
#'     \item{`g_wp`}{Combined conductance of wall and plasmalemma (rwp = rwall + rplasmalemma) in Mesophyll conductance in mol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`g_ch`}{Combined conductance of chloroplast envelope and stroma (rch = renvelope + rstroma) in Mesophyll conductance in mol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'   }
#' @return List of calculated rate:
#'   \describe{
#'     \item{`An`}{Net assimilation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Ac`}{Net assimilation, RuBP-saturated, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Aj`}{Net assimilation, RuBP regeneration-limited, in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`Limitation`}{Character vector specifying the limitation at which CO2 level.}
#'   }
#' @export


Tholen2012 = function(envs,pars) {

  calculateVc <- function(Cc, O, PAR,pars){
    pars$K_CO <- pars$K_C * (1 + O / pars$K_O)
    Wc <- pars$V_cmax * (Cc) / (Cc + pars$K_CO)
    bJ <- - (pars$alpha_J * PAR + pars$J_max)
    J <- (- bJ - sqrt((bJ)^2 - 4 * pars$theta_J * pars$alpha_J * PAR * pars$J_max)) / (2 * pars$theta_J)

    Wj <- (Cc * J) / (4 * Cc + 8 * pars$gamma_star)

    Vc <- min(Wc,Wj)
    limiting <- c("Wc", "Wj")[which.min(c(Wc, Wj))]
    return(list(Wc = Wc, Wj = Wj, Vc=Vc, lim_factor = limiting))
  }
  solve_for_Cc <- function(Ci, O, PAR, pars){
    fn_min <- function(Cc){
      pars$gamma_star <- O / (2 * pars$S_co)
      Vc <- calculateVc(Cc, O, PAR, pars)$Vc

      # Eq B1: Vo = Vc * gamma_star/Cc
      Vo <- Vc * pars$gamma_star/Cc

      # Eq 4 in the publication: A = Vc * (1 - gamma_star/Cc)- Rd
      An  <- Vc - Vo - pars$R_d

      # Eq 9 in the publication:
      gm <- pars$g_wp * pars$g_ch / (pars$g_wp + pars$g_ch * (1 + (Vo + pars$R_d)/An))

      # Eq 1 in the publication
      x <- Ci - Cc - An/gm
      return(x)
    }
    lo <- 1e-6; hi <- max(Ci, 1.2 * Ci + 5)
    # widen if needed
    f_lo <- fn_min(lo); f_hi <- fn_min(hi)
    if (f_lo * f_hi > 0) hi <- hi + 200
    Cc <- uniroot(fn_min, c(lo, hi), tol = 1e-9)$root
    return(Cc)
  }

  ret = list()
  ret$Limitation <- character(length(envs$C_i))
  ret$Wc <- numeric(length(envs$C_i))
  ret$Wj <- numeric(length(envs$C_i))
  ret$Wmin <- numeric(length(envs$C_i))
  ret$An <- numeric(length(envs$C_i))

  for (i in 1:length(envs$C_i)) {
    pars$gamma_star <- envs$O[i] / (2 * pars$S_co)
    Cc <- solve_for_Cc(envs$C_i[i],envs$O[i], envs$PPFD[i],pars)

    res <- calculateVc(envs$C_i[i],envs$O[i], envs$PPFD[i], pars)
    Vc <- res$Vc

    # Eq 4 in the publication: A = Vc * (1 - gamma_star/Cc)- Rd
    ret$An[i] <- Vc * (1 - pars$gamma_star/Cc)- pars$R_d

    ret$Wc[i] <- res$Wc
    ret$Wj[i] <- res$Wj
    ret$Wmin[i] <- res$Vc
    ret$Limitation <- res$lim_factor
  }

  gamma_star <- envs$O / (2 * pars$S_co)
  ret$Ac <- ret$Wc * (1 - gamma_star/Cc)- pars$R_d
  ret$Aj <- ret$Wj * (1 - gamma_star/Cc)- pars$R_d

  return(ret)
}

