#' Calculate carboxylation rate limited by three process
#' Wc: RuBP-saturated carboxylation rate
#' Wj: RuBP regeneration-limited carboxylation rate
#' Wp: TPU-limited carboxylation rate
#' @param envs List of environmental variables:
#'   \describe{
#'     \item{`PPFD`}{Photosynthetic photon flux density in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`C_i`}{CO2 concentration in chloroplast in μmol mol\eqn{^{-1}}.}
#'   }
#' @param pars List of parameters:
#'   \describe{
#'     \item{`V_cmax`}{Maximum rate of RuBP-saturated carboxylation in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`J_max`}{Maximum rate of electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`V_tpu`}{Rate of triose phosphate export from the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_C`}{Michaelis-Menten constant of Rubisco for CO2, same units as CO2.}
#'     \item{`K_O`}{Michaelis-Menten constant of Rubisco for O2, same units as O2.}
#'     \item{`S_co`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`max_alpha_G`}{Upper bound of glycolate carbon proportion taken out of the photorespiratory pathway as glycine.}
#'     \item{`max_alpha_S`}{Upper bound of glycolate carbon proportion taken out of the photorespiratory pathway as serine.}
#'     \item{`N_max`}{Maximum rate of de novo nitrogen supply to the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`alpha_J`}{quantum efficiency, unitless.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
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

  npoint <- length(envs$C_i)
  ret = list()
  ret$Limitation <- character(npoint)
  ret$Wc <- numeric(npoint)
  ret$Wj <- numeric(npoint)
  ret$Wp <- numeric(npoint)
  ret$Wmin <- numeric(npoint)

  ret$An <- numeric(npoint)

  gamma_star_agt <- matrix(nrow = 3, ncol = npoint)

  alpha_G_c <- numeric(npoint)
  alpha_G_j <- numeric(npoint)
  alpha_G_p <- numeric(npoint)


  alpha_G <- matrix( nrow = 3, ncol = npoint)
  alpha_S <- matrix( nrow = 3, ncol = npoint)

  ret$alpha_G <- numeric(npoint)
  ret$alpha_S <- numeric(npoint)
  for (i in 1:npoint) {
    Ci <- envs$C_i[i]
    PAR <- envs$PPFD[i]
    O <- envs$O[i]

    bJ <- - (pars$alpha_J * PAR + pars$J_max)
    J <- (- bJ - sqrt((bJ)^2 - 4 * pars$theta_J * pars$alpha_J * PAR * pars$J_max)) / (2 * pars$theta_J)
    # S_co <- (pars$V_cmax * pars$K_O) / (pars$V_omax * pars$K_C)
    # Beta: proportion of N being assimilated as glycine
    beta <- 3 * pars$max_alpha_G / (3*pars$max_alpha_G + 2*pars$max_alpha_S)

    phi <- O / ( pars$S_co * Ci)
    pars$K_CO <- pars$K_C * (1 + O / pars$K_O)

    # Calculate Vo, gamma_G and gamma_S under Rubisco limitation
    Vo_c <- pars$V_cmax * O / pars$S_co / (Ci + pars$K_CO)
    alpha_G[1,i] <- min(pars$max_alpha_G, pars$N_max * beta / Vo_c)
    alpha_S[1,i] <- min(pars$max_alpha_S, 3*pars$N_max * (1 - beta) / (2*Vo_c))

    if (J > (pars$N_max * (2*beta + 6))){
      alpha_G[2,i] <- min(pars$max_alpha_G,
                        4 * pars$N_max * beta * (1/phi + 1) / (J - pars$N_max * (2*beta + 6)))

      alpha_S[2,i] <- min(pars$max_alpha_S,
                        6*pars$N_max * (1 - beta) * (1/phi + 1) / (J - pars$N_max * (2*beta + 6)))

    }else{
      alpha_G[2,i] <- pars$max_alpha_G
      alpha_S[2,i] <- pars$max_alpha_S
    }

    # Calculate Vo, gamma_G and gamma_S under RuBP-regeneration limitation

    # Vo_j <- J / (4/phi + 4 + 8*alpha_G_j + 4*alpha_S)

    # Calculate gamma_G and gamma_S under TPU limitation

    alpha_G[3,i] <- min(pars$max_alpha_G, pars$N_max * beta * (2/phi - 1)/ (6*pars$V_tpu + 3*pars$N_max*(2-beta) ))
    alpha_S[3,i] <- min(pars$max_alpha_S, 3/2 * pars$N_max * (1 - beta) * (2/phi - 1) / (6*pars$V_tpu + 3*pars$N_max*(2-beta) ))

    ret$Wc[i] = pars$V_cmax * Ci / (Ci + pars$K_CO)
    ret$Wj[i] = J / (4 + (4 + 8*alpha_G[2,i] + 4*alpha_S[2,i]) * phi)
    ret$Wp[i] = 3 * pars$V_tpu / (1 - 0.5 * (1 + 3*alpha_G[3,i] + 4*alpha_S[3,i]) * phi)

    gamma_star_agt[1,i] <- 0.5 * (1-alpha_G[1,i]) * O / pars$S_co
    gamma_star_agt[2,i] <- 0.5 * (1-alpha_G[2,i]) * O / pars$S_co
    gamma_star_agt[3,i] <- 0.5 * (1-alpha_G[3,i]) * O / pars$S_co

    ret$Wmin[i] <- min(ret$Wc[i], ret$Wj[i], ret$Wp[i])
    min_ind <- which.min(c(ret$Wc[i], ret$Wj[i], ret$Wp[i]))

    ret$Limitation[i] <- c("Wc", "Wj", "Wp")[min_ind]
    ret$An[i] = (1- gamma_star_agt[min_ind,i]/Ci)*ret$Wmin[i] - pars$R_d
    ret$alpha_G[i] <- alpha_G[min_ind,i]
    ret$alpha_S[i] <- alpha_S[min_ind,i]


  }
  ret$Ac <- (1- gamma_star_agt[1,]/envs$C_i)*ret$Wc - pars$R_d
  ret$Aj <- (1- gamma_star_agt[2,]/envs$C_i)*ret$Wj - pars$R_d
  ret$Ap <- (1- gamma_star_agt[3,]/envs$C_i)*ret$Wp - pars$R_d



  return(ret)
}
