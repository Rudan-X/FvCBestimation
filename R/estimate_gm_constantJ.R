estimate_gm_constantJ <- function(A, Ci, pars){
  ACi <- data.frame(An = A, Ci = Ci)
  ACi <- ACi[ACi$Ci>500,]
  calc_varJ <- function(gm, ACi, pars) {
    Cc <- ACi$Ci - ACi$An / gm
    J <- (ACi$An + pars$R_d) * 4 * (Cc + 2 * pars$gamma_star) / (Cc - pars$gamma_star)

    var_J <- sum((mean(J) - J)^2 / (length(J) - 1))
    return(var_J)
  }

  # Optimize gm to minimize var_J
  opt_res <- optimize(
    f = calc_varJ,
    interval = c(0.01, 2),
    ACi = ACi,
    pars = pars
  )

  best_gm <- opt_res$minimum

  return(best_gm)
}
