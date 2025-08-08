#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' Wp: TPU-limited carboxylation rate
#' @param envs List of environmental variables:
#'   \describe{
#'     \item{`PPFD`}{Photosynthetic photon flux density in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`C_chl`}{CO2 concentration in chloroplast in μmol mol\eqn{^{-1}}.}
#'   }
#' @param pars List of parameters:
#'   \describe{
#'     \item{`V_cmax`}{Maximum rate of RuBP-saturated carboxylation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`J`}{Electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`V_tpu`}{Rate of triose phosphate export from the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as O2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`alpha_tpu`}{Fraction of  the glycolate carbon that is not returned to the chloroplast, unitless.}
#'     \item{`max_alpha_G`}{Upper bound of glycolate carbon proportion taken out of the photorespiratory pathway as glycine.}
#'     \item{`max_alpha_S`}{Upper bound of glycolate carbon proportion taken out of the photorespiratory pathway as serine.}
#'     \item{`N_max`}{Maximum rate of de novo nitrogen supply to the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'
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


Busch2018 = function(envs,pars) {

  S_co <- (pars$V_cmax * pars$K_O) / (pars$V_omax * pars$K_C)
  # Beta: proportion of N being assimilated as glycine
  beta <- 3 * pars$max_alpha_G / (3*pars$max_alpha_G + 2*pars$max_alpha_S)

  phi <- envs$O / ( S_co * envs$C_chl)
  pars$K_CO <- pars$K_C * (1 + envs$O / pars$K_O)

  # Calculate Vo, gamma_G and gamma_S under Rubisco limitation
  Vo_c <- pars$V_cmax * envs$O / S_co / (envs$C_chl + pars$K_CO)
  alpha_G_c <- pmin(pars$max_alpha_G, pars$N_max * beta / Vo_c)
  alpha_S_c <- pmin(pars$max_alpha_S, 3*pars$N_max * (1 - beta) / (2*Vo_c))

  # Calculate Vo, gamma_G and gamma_S under RuBP-regeneration limitation

  # Vo_j <- pars$J / (4/phi + 4 + 8*alpha_G_j + 4*alpha_S)

  if (pars$J > (pars$N_max * (2*beta + 6))){
    alpha_G_j <- pmin(pars$max_alpha_G,
                      4 * pars$N_max * beta * (1/phi + 1) / (pars$J - pars$N_max * (2*beta + 6)))

    alpha_S_j <- pmin(pars$max_alpha_S,
                      6*pars$N_max * (1 - beta) * (1/phi + 1) / (pars$J - pars$N_max * (2*beta + 6)))

  }else{
    alpha_G_j <- pars$max_alpha_G
    alpha_S_j <- pars$max_alpha_S
  }

  # Calculate gamma_G and gamma_S under TPU limitation

  alpha_G_p <- pmin(pars$max_alpha_G, pars$N_max * beta * (2/phi - 1)/ (6*pars$V_tpu + 3*pars$N_max*(2-beta) ))
  alpha_S_p <- pmin(pars$max_alpha_S, 3/2 * pars$N_max * (1 - beta) * (2/phi - 1) / (6*pars$V_tpu + 3*pars$N_max*(2-beta) ))

  ret = list(
    # Wc = pars$V_cmax * envs$C_chl / (envs$C_chl + pars$K_C * (1 + 1e3 * envs$O / pars$K_O)),
    Wc = pars$V_cmax * envs$C_chl / (envs$C_chl + pars$K_CO),
    Wj = pars$J / (4 + (4 + 8*alpha_G_j + 4*alpha_S_j) * phi),
    Wp = 3 * pars$V_tpu / (1 - 0.5 * (1 + 3*alpha_G_p + 4*alpha_S_p) * phi)
  )

  gamma_star_agt <- matrix(nrow = 3, ncol = length(envs$C_chl))

  gamma_star_agt[1,] <- 0.5 * (1-alpha_G_c) * envs$O / S_co
  gamma_star_agt[2,] <- 0.5 * (1-alpha_G_j) * envs$O / S_co
  gamma_star_agt[3,] <- 0.5 * (1-alpha_G_p) * envs$O / S_co


  alpha_G <- matrix(c(alpha_G_c, alpha_G_j,  alpha_G_p), nrow = 3, ncol = length(envs$C_chl), byrow = TRUE)
  alpha_S <- matrix(c(alpha_S_c, alpha_S_j,  alpha_S_p), nrow = 3, ncol = length(envs$C_chl), byrow = TRUE)


  ret$Wmin <- c()
  ret$An <- c()
  ret$alpha_G <- c()
  ret$alpha_S <- c()

  limiting_factor <- character(length(envs$C_chl))
  for (i in seq_along(envs$C_chl)){
    ret$Wmin[i] <- min(ret$Wc[i], ret$Wj[i], ret$Wp[i])
    min_ind <- which.min(c(ret$Wc[i], ret$Wj[i], ret$Wp[i]))

    limiting_factor[i] <- c("Wc", "Wj", "Wp")[min_ind]
    ret$An[i] = (1- gamma_star_agt[min_ind,i]/envs$C_chl[i])*ret$Wmin[i] - pars$R_d
    ret$alpha_G[i] <- alpha_G[min_ind,i]
    ret$alpha_S[i] <- alpha_S[min_ind,i]
  }

  ret$Ac <- (1- gamma_star_agt[1,]/envs$C_chl)*ret$Wc - pars$R_d
  ret$Aj <- (1- gamma_star_agt[2,]/envs$C_chl)*ret$Wj - pars$R_d
  ret$Ap <- (1- gamma_star_agt[3,]/envs$C_chl)*ret$Wp - pars$R_d


  ret$Limitation <- limiting_factor

  return(ret)
}
