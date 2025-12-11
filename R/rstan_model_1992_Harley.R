stan_Harley92 <- "
functions {
  real J_rh(real PPFD, real alpha_J, real J_max) {
    real J  = alpha_J * PPFD / sqrt(1 + square(alpha_J * PPFD / J_max));
    return J;
  }

  real Wc_from(real V_cmax, real C, real Kco) {
    return V_cmax * C / (C + Kco);
  }
  real Wj_from(real J, real C, real gamma_star) {
    return J * C / (4*C + 8*gamma_star);
  }

  // Solve Cc as quadratic root
  real solve_quad_physical(real a, real b, real c, real Ci) {
    real disc = b*b - 4*a*c;
    if (disc < 0) return negative_infinity(); // will be filtered later
    real sqrtD = sqrt(disc);
    // Stable formula
    real q = (b > 0) ? -0.5 * (b + sqrtD) : -0.5 * (b - sqrtD);
    real r1 = q / a;
    real r2 = c / q;
    // prefer a root in (0, Ci]
    real best = negative_infinity();
    if (r1 > 0 && r1 <= Ci) best = r1;
    if (r2 > 0 && r2 <= Ci) {
      if (best == negative_infinity()) best = r2;
      else if (fabs(r2 - Ci) < fabs(best - Ci)) best = r2; // heuristic
    }
    return best; // may be -inf if none valid
  }

  real model_Harley92(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real alpha_J, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;
    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real J = J_rh(PPFD, alpha_J, J_max);

    // Coefficients for Ac-branch
    real a_c = gm;
    real b_c = V_cmax + gm*(K_CO - Ci) - R_d;
    real c_c = -(K_CO*gm*Ci + K_CO*R_d + V_cmax*gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4*gm;
    real b_j = J + 8*gm*gamma_star - 4*gm*Ci - 4*R_d;
    real c_j = -(8*gm*gamma_star*Ci + 8*gamma_star*R_d + J*gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Ac = (1 - gamma_star / Cc_c) * Wc_from(V_cmax, Cc_c, K_CO) - R_d;
    real Aj = (1 - gamma_star / Cc_j) * Wj_from(J, Cc_j, gamma_star) - R_d;

    real An;

    if (Ac < Aj) {
      An  = Ac;
    } else {
      An  = Aj;
    }

    return An;
    }
}

data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // constants (estimated or fixed; you currently estimate them)
  // real K_C; real K_O; real alpha_J; real theta_J;
}
parameters {
  real<lower=50,  upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5>   R_d;
  real<lower=10,  upper=60>  gamma_star;
  real<lower=0.2, upper=0.5> alpha_J;
  // real<lower=0.3, upper=0.9> theta_J;
  real<lower=200, upper=300> K_C;
  real<lower=150, upper=300> K_O;

  real<lower=0.1, upper=0.7> gm;


  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Harley92(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_C, K_O, gamma_star, alpha_J,
                       gm);

  }

  // Likelihood
  A_obs ~ normal(A_hat, sigma);

  // Priors (examples)
  // V_cmax      ~ normal(125, 30);
  // J_max       ~ normal(175, 40);
  // R_d         ~ normal(2.5, 1.0);
  // gamma_star  ~ normal(35, 10);
  sigma       ~ normal(0, 2) T[0,];
}
"


stan_Harley92_B <- "
functions {
  real J_rh(real PPFD, real alpha_J, real J_max) {
    real J  = alpha_J * PPFD / sqrt(1 + square(alpha_J * PPFD / J_max));
    return J;
  }

  real Wc_from(real V_cmax, real C, real Kco) {
    return V_cmax * C / (C + Kco);
  }
  real Wj_from(real J, real C, real gamma_star) {
    return J * C / (4*C + 8*gamma_star);
  }

  // Solve Cc as quadratic root
  real solve_quad_physical(real a, real b, real c, real Ci) {
    real disc = b*b - 4*a*c;
    if (disc < 0) return negative_infinity(); // will be filtered later
    real sqrtD = sqrt(disc);
    // Stable formula
    real q = (b > 0) ? -0.5 * (b + sqrtD) : -0.5 * (b - sqrtD);
    real r1 = q / a;
    real r2 = c / q;
    // prefer a root in (0, Ci]
    real best = negative_infinity();
    if (r1 > 0 && r1 <= Ci) best = r1;
    if (r2 > 0 && r2 <= Ci) {
      if (best == negative_infinity()) best = r2;
      else if (fabs(r2 - Ci) < fabs(best - Ci)) best = r2; // heuristic
    }
    return best; // may be -inf if none valid
  }

  real model_Harley92(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real alpha_J, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;
    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real J = J_rh(PPFD, alpha_J, J_max);

    // Coefficients for Ac-branch
    real a_c = gm;
    real b_c = V_cmax + gm*(K_CO - Ci) - R_d;
    real c_c = -(K_CO*gm*Ci + K_CO*R_d + V_cmax*gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4*gm;
    real b_j = J + 8*gm*gamma_star - 4*gm*Ci - 4*R_d;
    real c_j = -(8*gm*gamma_star*Ci + 8*gamma_star*R_d + J*gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Ac = (1 - gamma_star / Cc_c) * Wc_from(V_cmax, Cc_c, K_CO) - R_d;
    real Aj = (1 - gamma_star / Cc_j) * Wj_from(J, Cc_j, gamma_star) - R_d;

    real An;

    if (Ac < Aj) {
      An  = Ac;
    } else {
      An  = Aj;
    }

    return An;
    }
}

data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // constants (estimated or fixed; you currently estimate them)
  real K_C;
  real K_O;
  real alpha_J;
  // real theta_J;
  real gm;
}
parameters {
  real<lower=50,  upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5>   R_d;
  real<lower=10,  upper=60>  gamma_star;

  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Harley92(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_C, K_O, gamma_star, alpha_J,
                       gm);

  }

  // Likelihood
  A_obs ~ normal(A_hat, sigma);

  // Priors (examples)
  // V_cmax      ~ normal(125, 30);
  // J_max       ~ normal(175, 40);
  // R_d         ~ normal(2.5, 1.0);
  // gamma_star  ~ normal(35, 10);
  sigma       ~ normal(0, 2) T[0,];
}
"

stan_Harley92_C <- "
functions {
  real J_rh(real PPFD, real alpha_J, real J_max) {
    real J  = alpha_J * PPFD / sqrt(1 + square(alpha_J * PPFD / J_max));
    return J;
  }

  real Wc_from(real V_cmax, real C, real Kco) {
    return V_cmax * C / (C + Kco);
  }
  real Wj_from(real J, real C, real gamma_star) {
    return J * C / (4*C + 8*gamma_star);
  }

  // Solve Cc as quadratic root
  real solve_quad_physical(real a, real b, real c, real Ci) {
    real disc = b*b - 4*a*c;
    if (disc < 0) return negative_infinity(); // will be filtered later
    real sqrtD = sqrt(disc);
    // Stable formula
    real q = (b > 0) ? -0.5 * (b + sqrtD) : -0.5 * (b - sqrtD);
    real r1 = q / a;
    real r2 = c / q;
    // prefer a root in (0, Ci]
    real best = negative_infinity();
    if (r1 > 0 && r1 <= Ci) best = r1;
    if (r2 > 0 && r2 <= Ci) {
      if (best == negative_infinity()) best = r2;
      else if (fabs(r2 - Ci) < fabs(best - Ci)) best = r2; // heuristic
    }
    return best; // may be -inf if none valid
  }

  real model_Harley92(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)

    real J = J_rh(PPFD, alpha_J, J_max);

    // Coefficients for Ac-branch
    real a_c = gm;
    real b_c = V_cmax + gm*(K_CO - Ci) - R_d;
    real c_c = -(K_CO*gm*Ci + K_CO*R_d + V_cmax*gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4*gm;
    real b_j = J + 8*gm*gamma_star - 4*gm*Ci - 4*R_d;
    real c_j = -(8*gm*gamma_star*Ci + 8*gamma_star*R_d + J*gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Ac = (1 - gamma_star / Cc_c) * Wc_from(V_cmax, Cc_c, K_CO) - R_d;
    real Aj = (1 - gamma_star / Cc_j) * Wj_from(J, Cc_j, gamma_star) - R_d;

    real An;

    if (Ac < Aj) {
      An  = Ac;
    } else {
      An  = Aj;
    }

    return An;
    }
}

data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // constants (estimated or fixed; you currently estimate them)
  // real K_C; real K_O; real alpha_J; real theta_J;
}
parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5>   R_d;
  real<lower=10,  upper=60>  gamma_star;
  real<lower=0.2, upper=0.5> alpha_J;
  // real<lower=0.3, upper=0.9> theta_J;
  real<lower=340, upper=720> K_CO;

  real<lower=0.1, upper=0.7> gm;


  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Harley92(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_CO, gamma_star, alpha_J,
                       gm);

  }

  // Likelihood
  A_obs ~ normal(A_hat, sigma);

  // Priors (examples)
  // V_cmax      ~ normal(125, 30);
  // J_max       ~ normal(175, 40);
  // R_d         ~ normal(2.5, 1.0);
  // gamma_star  ~ normal(35, 10);
  sigma       ~ normal(0, 2) T[0,];
}
"


stan_Harley92_D <- "
functions {
  real J_rh(real PPFD, real alpha_J, real J_max) {
    real J  = alpha_J * PPFD / sqrt(1 + square(alpha_J * PPFD / J_max));
    return J;
  }

  real Wc_from(real V_cmax, real C, real Kco) {
    return V_cmax * C / (C + Kco);
  }
  real Wj_from(real J, real C, real gamma_star) {
    return J * C / (4*C + 8*gamma_star);
  }

  // Solve Cc as quadratic root
  real solve_quad_physical(real a, real b, real c, real Ci) {
    real disc = b*b - 4*a*c;
    if (disc < 0) return negative_infinity(); // will be filtered later
    real sqrtD = sqrt(disc);
    // Stable formula
    real q = (b > 0) ? -0.5 * (b + sqrtD) : -0.5 * (b - sqrtD);
    real r1 = q / a;
    real r2 = c / q;
    // prefer a root in (0, Ci]
    real best = negative_infinity();
    if (r1 > 0 && r1 <= Ci) best = r1;
    if (r2 > 0 && r2 <= Ci) {
      if (best == negative_infinity()) best = r2;
      else if (fabs(r2 - Ci) < fabs(best - Ci)) best = r2; // heuristic
    }
    return best; // may be -inf if none valid
  }

  real model_Harley92(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)

    real J = J_rh(PPFD, alpha_J, J_max);

    // Coefficients for Ac-branch
    real a_c = gm;
    real b_c = V_cmax + gm*(K_CO - Ci) - R_d;
    real c_c = -(K_CO*gm*Ci + K_CO*R_d + V_cmax*gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4*gm;
    real b_j = J + 8*gm*gamma_star - 4*gm*Ci - 4*R_d;
    real c_j = -(8*gm*gamma_star*Ci + 8*gamma_star*R_d + J*gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Ac = (1 - gamma_star / Cc_c) * Wc_from(V_cmax, Cc_c, K_CO) - R_d;
    real Aj = (1 - gamma_star / Cc_j) * Wj_from(J, Cc_j, gamma_star) - R_d;

    real An;

    if (Ac < Aj) {
      An  = Ac;
    } else {
      An  = Aj;
    }

    return An;
    }
}

data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // constants (estimated or fixed; you currently estimate them)
  // real K_C; real K_O; real alpha_J; real theta_J;
}
parameters {
  real<lower=50, upper=400> V_cmax;
  real<lower=50, upper=400> J_max;
  real<lower=0.5, upper=5>   R_d;
  real<lower=10,  upper=60>  gamma_star;
  real<lower=0.2, upper=0.5> alpha_J;
  // real<lower=0.3, upper=0.9> theta_J;
  real<lower=300, upper=900> K_CO;

  real<lower=0.1, upper=0.7> gm;


  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Harley92(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_CO, gamma_star, alpha_J,
                       gm);

  }

  // Likelihood
  A_obs ~ normal(A_hat, sigma);

  // Priors (examples)
  // V_cmax      ~ normal(125, 30);
  // J_max       ~ normal(175, 40);
  // R_d         ~ normal(2.5, 1.0);
  // gamma_star  ~ normal(35, 10);
  sigma       ~ normal(0, 2) T[0,];
}

generated quantities {
  vector[N] log_lik;        // pointwise log-likelihood
  real log_lik_sum;         // total log-likelihood (convenience)

  for (n in 1:N) {
    real A_hat_n;
    A_hat_n = model_Harley92(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_CO, gamma_star, alpha_J,
                       gm);

    // pointwise log-likelihood for Normal(A_hat, sigma)
    log_lik[n] = normal_lpdf(A_obs[n] | A_hat_n, sigma);
  }

  log_lik_sum = sum(log_lik);
}
"
