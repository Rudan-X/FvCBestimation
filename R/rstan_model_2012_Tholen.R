stan_Tholen12 <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
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

  real model_Tholen12(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real alpha_J, real theta_J, real g_ch, real g_wp) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;
    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    // Coefficients for Ac-branch
    real a_c = g_ch;
    real b_c = - g_ch * Ci - R_d + g_ch * K_CO + V_cmax * (1 + g_ch / g_wp);
    real c_c = -(g_ch * Ci * K_CO + R_d * K_CO + V_cmax * gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4 * g_ch;
    real b_j = J * (1 + g_ch / g_wp) + 8 * g_ch * gamma_star - 4 * g_ch * Ci - 4 * R_d;
    real c_j = -(8 * g_ch * gamma_star * Ci + 8 * gamma_star * R_d + J * gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Wc = Wc_from(V_cmax, Cc_c, K_CO);
    real Wj = Wj_from(J, Cc_j, gamma_star);

    real Wmin;
    real Cstar;

    if (Wc < Wj) {
      Wmin  = Wc;
      Cstar = Cc_c;
    } else {
      Wmin  = Wj;
      Cstar = Cc_j;
    }

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=200, upper=300> K_C;
  real<lower=150, upper=300> K_O;
  real<lower=0.05, upper=0.7> g_ch;
  real<lower=0.1, upper=2> g_wp;


  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Tholen12(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_C, K_O, gamma_star, alpha_J, theta_J,
                       g_ch, g_wp);

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


stan_Tholen12_B <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
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

  real model_Tholen12(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real alpha_J, real theta_J, real g_ch, real g_wp) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;
    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    // Coefficients for Ac-branch
    real a_c = g_ch;
    real b_c = - g_ch * Ci - R_d + g_ch * K_CO + V_cmax * (1 + g_ch / g_wp);
    real c_c = -(g_ch * Ci * K_CO + R_d * K_CO + V_cmax * gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4 * g_ch;
    real b_j = J * (1 + g_ch / g_wp) + 8 * g_ch * gamma_star - 4 * g_ch * Ci - 4 * R_d;
    real c_j = -(8 * g_ch * gamma_star * Ci + 8 * gamma_star * R_d + J * gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Wc = Wc_from(V_cmax, Cc_c, K_CO);
    real Wj = Wj_from(J, Cc_j, gamma_star);

    real Wmin;
    real Cstar;

    if (Wc < Wj) {
      Wmin  = Wc;
      Cstar = Cc_c;
    } else {
      Wmin  = Wj;
      Cstar = Cc_j;
    }

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  real theta_J;
  real g_ch;
  real g_wp;

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

    A_hat[n]  = model_Tholen12(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_C, K_O, gamma_star, alpha_J, theta_J,
                       g_ch, g_wp);

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


stan_Tholen12_C <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
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

  real model_Tholen12(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J, real g_ch, real g_wp) {

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    // Coefficients for Ac-branch
    real a_c = g_ch;
    real b_c = - g_ch * Ci - R_d + g_ch * K_CO + V_cmax * (1 + g_ch / g_wp);
    real c_c = -(g_ch * Ci * K_CO + R_d * K_CO + V_cmax * gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4 * g_ch;
    real b_j = J * (1 + g_ch / g_wp) + 8 * g_ch * gamma_star - 4 * g_ch * Ci - 4 * R_d;
    real c_j = -(8 * g_ch * gamma_star * Ci + 8 * gamma_star * R_d + J * gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Wc = Wc_from(V_cmax, Cc_c, K_CO);
    real Wj = Wj_from(J, Cc_j, gamma_star);

    real Wmin;
    real Cstar;

    if (Wc < Wj) {
      Wmin  = Wc;
      Cstar = Cc_c;
    } else {
      Wmin  = Wj;
      Cstar = Cc_j;
    }

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=340, upper=720> K_CO;
  real<lower=0.05, upper=0.7> g_ch;
  real<lower=0.1, upper=2> g_wp;


  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Tholen12(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_CO, gamma_star, alpha_J, theta_J,
                       g_ch, g_wp);

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

stan_Tholen12_D <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
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

  real model_Tholen12(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J, real g_ch, real g_wp) {

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    // Coefficients for Ac-branch
    real a_c = g_ch;
    real b_c = - g_ch * Ci - R_d + g_ch * K_CO + V_cmax * (1 + g_ch / g_wp);
    real c_c = -(g_ch * Ci * K_CO + R_d * K_CO + V_cmax * gamma_star);

    // Coefficients for Aj-branch
    real a_j = 4 * g_ch;
    real b_j = J * (1 + g_ch / g_wp) + 8 * g_ch * gamma_star - 4 * g_ch * Ci - 4 * R_d;
    real c_j = -(8 * g_ch * gamma_star * Ci + 8 * gamma_star * R_d + J * gamma_star);

    real Cc_c = solve_quad_physical(a_c, b_c, c_c, Ci);
    real Cc_j = solve_quad_physical(a_j, b_j, c_j, Ci);


    real Wc = Wc_from(V_cmax, Cc_c, K_CO);
    real Wj = Wj_from(J, Cc_j, gamma_star);

    real Wmin;
    real Cstar;

    if (Wc < Wj) {
      Wmin  = Wc;
      Cstar = Cc_c;
    } else {
      Wmin  = Wj;
      Cstar = Cc_j;
    }

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=300, upper=900> K_CO;
  real<lower=0.05, upper=0.7> g_ch;
  real<lower=0.1, upper=2> g_wp;


  real<lower=0> sigma;
}
model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {

    A_hat[n]  = model_Tholen12(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_CO, gamma_star, alpha_J, theta_J,
                       g_ch, g_wp);

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
    A_hat_n = model_Tholen12(Ci[n], O[n], PPFD[n],
                       V_cmax, J_max, R_d,
                       K_CO, gamma_star, alpha_J, theta_J,
                       g_ch, g_wp);

    // pointwise log-likelihood for Normal(A_hat, sigma)
    log_lik[n] = normal_lpdf(A_obs[n] | A_hat_n, sigma);
  }

  log_lik_sum = sum(log_lik);
}
"
