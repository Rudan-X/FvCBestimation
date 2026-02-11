
J_nrh <- function(PPFD, alpha_J, theta_J, J_max) {
  b <- -(alpha_J * PPFD + J_max)
  disc <- pmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max)
  (-b - sqrt(disc)) / (2 * theta_J)
}

J_rh <- function(PPFD, alpha_J, J_max) {
  alpha_J * PPFD / sqrt(1 + (alpha_J * PPFD / J_max)^2)
}

Wc_from <- function(V_cmax, C, Kco) V_cmax * C / (C + Kco)

Wj_from <- function(J, C, gamma_star) J * C / (4*C + 8*gamma_star)

Wp_from <- function(V_tpu, C, gamma_star, alpha_G) {
  denom <- C - (1 + 3*alpha_G) * gamma_star
  if (denom <= 0) Inf else 3 * V_tpu * C / denom
}

Cref <- function(Cc, Ci, has_Cc) if (has_Cc == 1) Cc else Ci

# Solve Cc as quadratic root, selecting a "physical" root in (0, Ci)
solve_quad_physical <- function(a, b, c, Ci) {
  root <- c()
  for (x in 1:length(Ci)){
    disc <- b^2 - 4*a*c
    if (disc < 0) return(NA_real_)  # no real roots

    sqrtD <- sqrt(disc)
    # numerically stable quadratic formula
    if (b > 0) {
      q <- -0.5 * (b + sqrtD)
    } else {
      q <- -0.5 * (b - sqrtD)
    }
    r1 <- q / a
    r2 <- c / q

    # keep only roots in (0, Ci]
    valid_roots <- c(r1, r2)
    valid_roots <- valid_roots[valid_roots > 0 & valid_roots <= Ci]

    if (length(valid_roots) == 0) return(NA_real_)

    # heuristic: choose the root closest to Ci
    root <- valid_roots[which.min(abs(valid_roots - Ci))]
  }

  return(root)
}

solve_quad_physical_vecCi <- function(a, b, c, Ci) {
  mapply(function(ai, bi, ci, Cii) solve_quad_physical(ai, bi, ci, Cii),
         a, b, c, Ci)
}

model_FvCB80 <- function(Ci, has_Cc, Cc, O, PPFD,
                            V_cmax, J_max, R_d,
                            K_C, K_O, gamma_star,
                            alpha_J, theta_J) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000

  Cstar <- Cref(Cc, Ci, has_Cc)
  K_CO  <- K_C * (1 + O_conv / Ko_conv)

  J    <- J_nrh(PPFD, alpha_J, theta_J, J_max)
  Wc   <- Wc_from(V_cmax, Cstar, K_CO)
  Wj   <- Wj_from(J, Cstar, gamma_star)
  Wmin <- pmin(Wc, Wj)
  lim <- ifelse(Wc <= Wj, "Wc", "Wj")

  An <- (1 - gamma_star / Cstar) * Wmin - R_d
  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)

}

model_Harley92 <- function(Ci, O, PPFD,
                              V_cmax, J_max, R_d,
                              K_C, K_O, gamma_star,
                              alpha_J, gm) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000
  K_CO    <- K_C * (1 + O_conv / Ko_conv)

  J <- J_rh(PPFD, alpha_J, J_max)

  # Ac branch quadratic coefficients
  a_c <- gm
  b_c <- V_cmax + gm*(K_CO - Ci) - R_d
  c_c <- -(K_CO*gm*Ci + K_CO*R_d + V_cmax*gamma_star)

  # Aj branch quadratic coefficients
  a_j <- 4*gm
  b_j <- J + 8*gm*gamma_star - 4*gm*Ci - 4*R_d
  c_j <- -(8*gm*gamma_star*Ci + 8*gamma_star*R_d + J*gamma_star)

  Cc_c <- solve_quad_physical_vecCi(a_c, b_c, c_c, Ci)
  Cc_j <- solve_quad_physical_vecCi(a_j, b_j, c_j, Ci)

  Ac <- (1 - gamma_star / Cc_c) * Wc_from(V_cmax, Cc_c, K_CO) - R_d
  Aj <- (1 - gamma_star / Cc_j) * Wj_from(J, Cc_j, gamma_star) - R_d

  An <- pmin(Ac, Aj)
  lim <- ifelse(Ac <= Aj, "Ac", "Aj")

  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)
}

model_Tholen12 <- function(Ci, O, PPFD,
                              V_cmax, J_max, R_d,
                              K_C, K_O, gamma_star,
                              alpha_J, theta_J, g_ch, g_wp) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000
  K_CO    <- K_C * (1 + O_conv / Ko_conv)

  J <- J_nrh(PPFD, alpha_J, theta_J, J_max)

  # Ac branch
  a_c <- g_ch
  b_c <- - g_ch * Ci - R_d + g_ch * K_CO + V_cmax * (1 + g_ch / g_wp)
  c_c <- -(g_ch * Ci * K_CO + R_d * K_CO + V_cmax * gamma_star)

  # Aj branch
  a_j <- 4 * g_ch
  b_j <- J * (1 + g_ch / g_wp) + 8 * g_ch * gamma_star - 4 * g_ch * Ci - 4 * R_d
  c_j <- -(8 * g_ch * gamma_star * Ci + 8 * gamma_star * R_d + J * gamma_star)
  Cc_c <- solve_quad_physical(a_c, b_c, c_c, Ci)
  Cc_j <- solve_quad_physical(a_j, b_j, c_j, Ci)


  Wc <- Wc_from(V_cmax, Cc_c, K_CO)
  Wj <- Wj_from(J, Cc_j, gamma_star)

  is_c <- Wc <= Wj

  Wmin  <- ifelse(is_c, Wc, Wj)
  Cstar <- ifelse(is_c, Cc_c, Cc_j)
  lim <- ifelse(is_c, "Wc", "Wj")

  An <- (1 - gamma_star / Cstar) * Wmin - R_d
  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)
}

model_Caemmerer00 <- function(Ci, has_Cc, Cc, O, PPFD,
                                 V_cmax, J_max, R_d,
                                 K_C, K_O, gamma_star,
                                 alpha_J, theta_J, V_tpu, alpha_G) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000

  Cstar <- Cref(Cc, Ci, has_Cc)
  K_CO  <- K_C * (1 + O_conv / Ko_conv)

  J <- J_nrh(PPFD, alpha_J, theta_J, J_max)

  Wc   <- Wc_from(V_cmax, Cstar, K_CO)
  Wj   <- Wj_from(J, Cstar, gamma_star)
  Wp   <- Wp_from(V_tpu, Cstar, gamma_star, alpha_G)
  Wmin <- pmin(Wc, pmin(Wj, Wp))
  lim <- ifelse(Wc <= Wj & Wc <= Wp, "Wc",
                ifelse(Wj <= Wp, "Wj", "Wp"))

  An <- (1 - gamma_star / Cstar) * Wmin - R_d
  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)
}

model_Ethier04 <- function(Ci, O, PPFD,
                              V_cmax, J_max, R_d,
                              K_C, K_O, gamma_star,
                              alpha_J, theta_J, gm) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000
  K_CO    <- K_C * (1 + O_conv / Ko_conv)

  J <- J_nrh(PPFD, alpha_J, theta_J, J_max)

  a1 <- -1 / gm
  b1 <- ((V_cmax - R_d) / gm) + Ci + K_CO
  c1 <- R_d * (Ci + K_CO) - V_cmax * (Ci - gamma_star)
  Ac <- (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1)

  a2 <- -4 / gm
  b2 <- 4 * (Ci + 2 * gamma_star) - 4 * R_d/gm + J/gm
  c2 <- 4 * R_d * (Ci + 2 * gamma_star) - J * (Ci - gamma_star)
  Aj <- (-b2 + sqrt(b2*b2 - 4*a2*c2)) / (2*a2)

  An <- pmin(Ac, Aj)
  lim <- ifelse(Ac <= Aj, "Ac", "Aj")

  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)

}

model_Dubois07 <- function(Ci, has_Cc, Cc, O, PPFD,
                              V_cmax, J_max, R_d,
                              K_C, K_O, gamma_star,
                              alpha_J, theta_J, V_tpu, alpha_G) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000

  Cstar <- Cref(Cc, Ci, has_Cc)
  K_CO  <- K_C * (1 + O_conv / Ko_conv)

  J <- J_nrh(PPFD, alpha_J, theta_J, J_max)

  Ac <- (1 - gamma_star / Cstar) * Wc_from(V_cmax, Cstar, K_CO) - R_d
  Aj <- (1 - gamma_star / Cstar) * Wj_from(J, Cstar, gamma_star) - R_d
  Ap <- (1 - gamma_star / Cstar) * Wp_from(V_tpu, Cstar, gamma_star, alpha_G) - R_d

  An <- pmin(Ac, pmin(Aj, Ap))
  lim <- ifelse(Ac <= Aj & Ac <= Ap, "Ac",
                ifelse(Aj <= Ap, "Aj", "Ap"))

  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)

}

model_Yin04 <- function(Ci, has_Cc, Cc, O, PPFD,
                           V_cmax, J_max, R_d,
                           K_C, K_O, gamma_star,
                           phi2m, theta_J, f_Q, f_pseudo, h) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000

  Cstar <- Cref(Cc, Ci, has_Cc)
  K_CO  <- K_C * (1 + O_conv / Ko_conv)

  Wc <- Wc_from(V_cmax, Cstar, K_CO)

  a <- 4*Cstar + 8*gamma_star
  b <- 3*Cstar + 7*gamma_star

  denom <- a/(h*b) - 1
  if (abs(denom) < 1e-12) denom <- if (denom >= 0) 1e-12 else -1e-12

  f_cyc <- ( a*(2 + f_Q)/(h*b) - (1 - f_pseudo) ) / denom
  f_cyc <- pmin(1, pmax(0, f_cyc))

  one_minus_fc <- pmax(0, 1 - f_cyc)
  alpha_J <- one_minus_fc / (1 + (one_minus_fc / pmax(phi2m, 1e-12)))

  J <- J_nrh(PPFD, alpha_J, theta_J, J_max)

  Wj <- J * ((2 + f_Q - f_cyc) * Cstar) /
    (h * (3*Cstar + 7*gamma_star) * pmax(1e-12, (1 - f_cyc)))

  Wmin <- pmin(Wc, Wj)
  lim <- ifelse(Wc <= Wj, "Wc", "Wj")

  An <- (1 - gamma_star / Cstar) * Wmin - R_d

  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)
}

# ---------- Busch 2018 alpha helpers ----------
alpha_beta <- function(max_aG, max_aS) 3*max_aG / (3*max_aG + 2*max_aS)

phi_from <- function(Cstar, gamma_star) 2 * gamma_star / pmax(Cstar, 1e-9)

alpha_G_c <- function(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS) {
  Vo_c <- Vcmax * (2*gamma_star) / (Cstar + K_CO)
  beta <- alpha_beta(max_aG, max_aS)
  val  <- N_max * beta / pmax(Vo_c, 1e-12)
  pmin(max_aG, pmax(0, val))
}

alpha_S_c <- function(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS) {
  Vo_c <- Vcmax * (2*gamma_star) / (Cstar + K_CO)
  beta <- alpha_beta(max_aG, max_aS)
  val  <- (3 * N_max * (1 - beta)) / (2 * pmax(Vo_c, 1e-12))
  pmin(max_aS, pmax(0, val))
}

alpha_G_j <- function(Cstar, J, gamma_star, N_max, max_aG, max_aS) {
  beta <- alpha_beta(max_aG, max_aS)
  phi  <- phi_from(Cstar, gamma_star)
  thresh <- N_max * (2*beta + 6)
  if (J <= thresh) return(max_aG)
  val <- 4 * N_max * beta * (1/phi + 1) / (J - thresh)
  pmin(max_aG, pmax(0, val))
}

alpha_S_j <- function(Cstar, J, gamma_star, N_max, max_aG, max_aS) {
  beta <- alpha_beta(max_aG, max_aS)
  phi  <- phi_from(Cstar, gamma_star)
  thresh <- N_max * (2*beta + 6)
  if (J <= thresh) return(max_aS)
  val <- 6 * N_max * (1 - beta) * (1/phi + 1) / (J - thresh)
  pmin(max_aS, pmax(0, val))
}

alpha_G_p <- function(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS) {
  beta <- alpha_beta(max_aG, max_aS)
  phi  <- phi_from(Cstar, gamma_star)
  num  <- N_max * beta * pmax(2/phi - 1, 0)
  den  <- 6 * Vtpu + 3 * N_max * (2 - beta)
  val  <- num / pmax(den, 1e-12)
  pmin(max_aG, pmax(0, val))
}

alpha_S_p <- function(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS) {
  beta <- alpha_beta(max_aG, max_aS)
  phi  <- phi_from(Cstar, gamma_star)
  num  <- 1.5 * N_max * (1 - beta) * pmax(2/phi - 1, 0)
  den  <- 6 * Vtpu + 3 * N_max * (2 - beta)
  val  <- num / pmax(den, 1e-12)
  pmin(max_aS, pmax(0, val))
}

model_Busch18 <- function(Ci, has_Cc, Cc,
                             O, PPFD,
                             Vcmax, J_max, R_d,
                             K_C, K_O, gamma_star,
                             alpha_J, theta_J, Vtpu,
                             N_max, max_aG, max_aS) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000

  K_CO  <- K_C * (1 + O_conv / Ko_conv)
  Cstar <- Cref(Cc, Ci, has_Cc)

  J   <- J_nrh(PPFD, alpha_J, theta_J, J_max)
  Wc  <- Wc_from(Vcmax, Cstar, K_CO)

  aGj <- alpha_G_j(Cstar, J, gamma_star, N_max, max_aG, max_aS)
  aSj <- alpha_S_j(Cstar, J, gamma_star, N_max, max_aG, max_aS)
  phi <- phi_from(Cstar, gamma_star)
  Wj  <- J / (4 + (4 + 8*aGj + 4*aSj) * phi)

  aGp <- alpha_G_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS)
  aSp <- alpha_S_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS)
  denom <- 1 - 0.5 * (1 + 3*aGp + 4*aSp) * phi
  Wp <- if (denom <= 0) Inf else 3 * Vtpu / denom

  Wmin <- pmin(Wc, pmin(Wj, Wp))
  lim <- ifelse(Wc <= Wj & Wc <= Wp, "Wc",
                ifelse(Wj <= Wp, "Wj", "Wp"))

  g_c <- (1 - alpha_G_c(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS)) * gamma_star
  g_j <- (1 - aGj) * gamma_star
  g_p <- (1 - aGp) * gamma_star

  # mimic Stan's exact equality comparisons (can be fragile numerically)
  g_pick <- if (Wmin == Wc) g_c else if (Wmin == Wj) g_j else g_p

  An <- (1 - g_pick / Cstar) * Wmin - R_d
  df <- list(
    An = An,
    Limitation = lim
  )
  return(df)
}

# Xiao20 returns a length-2 vector: c(An, Phi2)
model_Xiao20 <- function(Ci, O, PPFD,
                            V_cmax, J_max, R_d,
                            K_C, K_O, gamma_star,
                            theta_J, s, Phi2LL, gm) {
  O_conv  <- O * 1000
  Ko_conv <- K_O * 1000

  K_CO <- K_C * (1 + O_conv / Ko_conv)

  Ji <- PPFD * s * Phi2LL
  bJ <- -(Ji + J_max)
  J  <- (-bJ - sqrt(bJ^2 - 4 * theta_J * Ji * J_max)) /(2 * theta_J)

  a1 <- -1 / gm
  b1 <- ((V_cmax - R_d) / gm) + Ci + K_CO
  c1 <- R_d * (Ci + K_CO) - V_cmax * (Ci - gamma_star)
  Ac <- (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1)

  a2 <- -4 / gm
  b2 <- 4 * (Ci + 2 * gamma_star) - 4 * R_d/gm + J/gm
  c2 <- 4 * R_d * (Ci + 2 * gamma_star) - J * (Ci - gamma_star)
  Aj <- (-b2 + sqrt(b2*b2 - 4*a2*c2)) / (2*a2)

  An <- pmin(Ac, Aj)
  lim <- ifelse(Ac <= Aj, "Ac", "Aj")


  Cc <- Ci - An/gm
  J_a <- ((An + R_d) * (4*Cc + 8*gamma_star)) / (Cc - gamma_star)
  J_f <- pmin(J_a, J_max)
  Phi2 <- J_f / (PPFD * s)

  df <- list(
    An = An,
    Phi2= Phi2,
    Limitation = lim
  )
  return(df)
}
