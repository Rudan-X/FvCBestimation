#' Simulate FvCB1980
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
#'     \item{`J_max`}{Maximum rate of electron transport through PSII in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`V_tpu`}{Rate of triose phosphate export from the chloroplast in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`K_CO`}{Apparent Michaelis-Menten constant of Rubisco, same units as CO2.}
#'     \item{`R_d`}{Day respiration in μmol m\eqn{^{-2}} s\eqn{^{-1}}.}
#'     \item{`gamma_star`}{CO2 photocompensation point, same units as CO2.}
#'     \item{`phi_J`}{quantum efficiency, unitless.}
#'     \item{`theta_J`}{Convexity factor for response of J to light, unitless.}
#'   }
#' @return List of two ggplot graphs:
#'   \describe{
#'     \item{`An`} Graph showing net assimilation (An) versus Ci
#'     \item{`Wmin`} Graph showing minimum between W_c, W_j, and W_p versus Ci
#'   }
#' @export
simulate_FvCB1980 <- function(envs,pars){
  Anet=FvCB1980(envs,pars)

  df <- as.data.frame(Anet)
  df$Ci <- envs$C_chl # add a time or observation index
  df <- df[,c("Ac","Aj","Ap","An","Ci")]

  # Reshape to long format
  df_long <- df %>%
    pivot_longer(cols = -Ci, names_to = "Limitation", values_to = "Value")


  df_long <- df_long %>%
    mutate(Linetype = ifelse(Limitation == "An", "solid", "dashed"))

  df_long$Limitation <- factor(df_long$Limitation, levels=c("Ac","Aj","Ap","An"))
  # Plot
  g1 <- ggplot(df_long, aes(x = Ci, y = Value, color = Limitation, linetype = Linetype)) +
    geom_line(size = 1, alpha = 0.9) +
    scale_linetype_identity() +  # Use values from the Linetype column directly
    theme_minimal(base_size = 14) +
    labs(
      x = expression(C[i]~"(µmol mol"^{-1}*")"),
      y = expression(An~"(µmol m"^{-2}~s^{-1}*")"),
      color = NULL
    ) +
    theme(
      legend.position = "right"
    )


  df <- as.data.frame(Anet)
  df$Ci <- envs$C_chl # add a time or observation index
  df <- df[,c(1,2,3,4,10)]

  # Reshape to long format
  df_long <- df %>%
    pivot_longer(cols = -Ci, names_to = "Limitation", values_to = "Value")


  df_long <- df_long %>%
    mutate(Linetype = ifelse(Limitation == "Wmin", "solid", "dashed"))

  df_long$Limitation <- factor(df_long$Limitation, levels=c("Wc","Wj","Wp","Wmin"))
  # Plot
  g2 <- ggplot(df_long, aes(x = Ci, y = Value, color = Limitation, linetype = Linetype)) +
    geom_line(size = 1, alpha = 0.9) +
    scale_linetype_identity() +  # Use values from the Linetype column directly
    theme_minimal(base_size = 14) +
    labs(
      x = expression(C[i]~"(µmol mol"^{-1}*")"),
      y = expression(W~"(µmol m"^{-2}~s^{-1}*")"),
      color = NULL
    ) +
    theme(
      legend.position = "right"
    )

  return(list(An=g1, Wmin=g2))
}
