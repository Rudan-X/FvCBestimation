stan_Xiao20_An <- "
functions {
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Xiao20(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real theta_J, real s, real Phi2LL, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real Ji = PPFD * s * Phi2LL;
    real bJ = -(Ji + J_max);
    real J = (-bJ - sqrt(bJ^2 - 4 * theta_J * Ji * J_max)) /(2 * theta_J);

    real a1 = -1.0 / gm;
    real b1 = ((V_cmax - R_d) / gm) + Ci + K_CO;
    real c1 = R_d * (Ci + K_CO) - V_cmax * (Ci - gamma_star);
    real Ac = (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1);

    real a2 = -4.0 / gm;
    real b2 = 4 * (Ci + 2 * gamma_star) - 4 * R_d/gm + J/gm;
    real c2 = 4 * R_d * (Ci + 2 * gamma_star) - J * (Ci - gamma_star);
    real Aj = (-b2 + sqrt(b2*b2 - 4*a2*c2)) / (2*a2);

    return fmin(Ac, Aj);
  }
}
data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // Fixed parameters (constants passed from R)
  // real K_C;
  // real K_O;
  // real theta_J;
}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;
  real<lower=0.1,upper=0.7>  gm;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=200, upper=300> K_C;
  real<lower=150, upper=300> K_O;
  real<lower=0.3,upper=0.6>  s;          // leaf absorptance/partitioning factor
  real<lower=0.5,upper=0.9> Phi2LL;     // lower-level Phi2 driver (e.g., instrument LL)

  real<lower=0> sigma_A;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
     A_hat[n] = model_Xiao20(Ci[n], O[n], PPFD[n],
                                    V_cmax, J_max, R_d,
                                    K_C, K_O, gamma_star,
                                    theta_J, s, Phi2LL, gm);
  }

  A_obs ~ normal(A_hat, sigma_A);
  sigma_A     ~ normal(0, 2) T[0,];

}


"

stan_Xiao20_An_B <- "
functions {
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Xiao20(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real theta_J, real s, real Phi2LL, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real Ji = PPFD * s * Phi2LL;
    real bJ = -(Ji + J_max);
    real J = (-bJ - sqrt(bJ^2 - 4 * theta_J * Ji * J_max)) /(2 * theta_J);

    real a1 = -1.0 / gm;
    real b1 = ((V_cmax - R_d) / gm) + Ci + K_CO;
    real c1 = R_d * (Ci + K_CO) - V_cmax * (Ci - gamma_star);
    real Ac = (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1);

    real a2 = -4.0 / gm;
    real b2 = 4 * (Ci + 2 * gamma_star) - 4 * R_d/gm + J/gm;
    real c2 = 4 * R_d * (Ci + 2 * gamma_star) - J * (Ci - gamma_star);
    real Aj = (-b2 + sqrt(b2*b2 - 4*a2*c2)) / (2*a2);

    return fmin(Ac, Aj);
  }
}
data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // Fixed parameters (constants passed from R)
  real K_C;
  real K_O;
  real theta_J;
  real gm;
  real s;
  real Phi2LL;
}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;


  real<lower=0> sigma_A;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
     A_hat[n] = model_Xiao20(Ci[n], O[n], PPFD[n],
                                    V_cmax, J_max, R_d,
                                    K_C, K_O, gamma_star,
                                    theta_J, s, Phi2LL, gm);
  }

  A_obs ~ normal(A_hat, sigma_A);
  sigma_A     ~ normal(0, 2) T[0,];

}
"

stan_Xiao20_An_C <- "
functions {
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Xiao20(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real theta_J, real s, real Phi2LL, real gm) {

    real Ji = PPFD * s * Phi2LL;
    real bJ = -(Ji + J_max);
    real J = (-bJ - sqrt(bJ^2 - 4 * theta_J * Ji * J_max)) /(2 * theta_J);

    real a1 = -1.0 / gm;
    real b1 = ((V_cmax - R_d) / gm) + Ci + K_CO;
    real c1 = R_d * (Ci + K_CO) - V_cmax * (Ci - gamma_star);
    real Ac = (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1);

    real a2 = -4.0 / gm;
    real b2 = 4 * (Ci + 2 * gamma_star) - 4 * R_d/gm + J/gm;
    real c2 = 4 * R_d * (Ci + 2 * gamma_star) - J * (Ci - gamma_star);
    real Aj = (-b2 + sqrt(b2*b2 - 4*a2*c2)) / (2*a2);

    return fmin(Ac, Aj);
  }
}
data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // Fixed parameters (constants passed from R)
  // real K_C;
  // real K_O;
  // real theta_J;
}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;
  real<lower=0.1,upper=0.7>  gm;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=340, upper=720> K_CO;
  real<lower=0.3,upper=0.6>  s;          // leaf absorptance/partitioning factor
  real<lower=0.5,upper=0.9> Phi2LL;     // lower-level Phi2 driver (e.g., instrument LL)

  real<lower=0> sigma_A;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
     A_hat[n] = model_Xiao20(Ci[n], O[n], PPFD[n],
                                    V_cmax, J_max, R_d,
                                    K_CO, gamma_star,
                                    theta_J, s, Phi2LL, gm);
  }

  A_obs ~ normal(A_hat, sigma_A);
  sigma_A     ~ normal(0, 2) T[0,];

}
"

stan_Xiao20_An_D <- "
functions {
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Xiao20(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real theta_J, real s, real Phi2LL, real gm) {

    real Ji = PPFD * s * Phi2LL;
    real bJ = -(Ji + J_max);
    real J = (-bJ - sqrt(bJ^2 - 4 * theta_J * Ji * J_max)) /(2 * theta_J);

    real a1 = -1.0 / gm;
    real b1 = ((V_cmax - R_d) / gm) + Ci + K_CO;
    real c1 = R_d * (Ci + K_CO) - V_cmax * (Ci - gamma_star);
    real Ac = (-b1 + sqrt(b1*b1 - 4*a1*c1)) / (2*a1);

    real a2 = -4.0 / gm;
    real b2 = 4 * (Ci + 2 * gamma_star) - 4 * R_d/gm + J/gm;
    real c2 = 4 * R_d * (Ci + 2 * gamma_star) - J * (Ci - gamma_star);
    real Aj = (-b2 + sqrt(b2*b2 - 4*a2*c2)) / (2*a2);

    return fmin(Ac, Aj);
  }
}
data {
  int<lower=1> N;
  vector[N] Ci;
  vector[N] O;
  vector[N] PPFD;
  vector[N] A_obs;

  // Fixed parameters (constants passed from R)
  // real K_C;
  // real K_O;
  // real theta_J;
}

parameters {
  real<lower=50, upper=400> V_cmax;
  real<lower=50, upper=400> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;
  real<lower=0.1,upper=0.7>  gm;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=300, upper=900> K_CO;
  real<lower=0.3,upper=0.6>  s;          // leaf absorptance/partitioning factor
  real<lower=0.5,upper=0.9> Phi2LL;     // lower-level Phi2 driver (e.g., instrument LL)

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
     A_hat[n] = model_Xiao20(Ci[n], O[n], PPFD[n],
                                    V_cmax, J_max, R_d,
                                    K_CO, gamma_star,
                                    theta_J, s, Phi2LL, gm);
  }

  A_obs ~ normal(A_hat, sigma);
  sigma     ~ normal(0, 2) T[0,];

}

generated quantities {
  vector[N] log_lik;        // pointwise log-likelihood
  real log_lik_sum;         // total log-likelihood (convenience)

  for (n in 1:N) {
    real A_hat_n;
    A_hat_n = model_Xiao20(Ci[n], O[n], PPFD[n],
                                    V_cmax, J_max, R_d,
                                    K_CO, gamma_star,
                                    theta_J, s, Phi2LL, gm);

    // pointwise log-likelihood for Normal(A_hat, sigma)
    log_lik[n] = normal_lpdf(A_obs[n] | A_hat_n, sigma);
  }

  log_lik_sum = sum(log_lik);
}
"
