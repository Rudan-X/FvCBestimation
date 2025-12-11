stan_Ethier04 <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real Ethier_model(real Ci,real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real alpha_J, real theta_J, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

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
  // real alpha_J;
  // real theta_J;

}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;
  real<lower=0.1,upper=0.7>  gm;
  real<lower=0.2, upper=0.5> alpha_J;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=200, upper=300> K_C;
  real<lower=150, upper=300> K_O;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = Ethier_model(Ci[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_C, K_O, gamma_star,
                          alpha_J, theta_J, gm);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // gm ~ normal(0.4, 0.15);
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}
"

stan_Ethier04_B <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real Ethier_model(real Ci,real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real alpha_J, real theta_J, real gm) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

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
  real alpha_J;
  real theta_J;
  real gm;
}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = Ethier_model(Ci[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_C, K_O, gamma_star,
                          alpha_J, theta_J, gm);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // gm ~ normal(0.4, 0.15);
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}
"

stan_Ethier04_C <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real Ethier_model(real Ci,real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J, real gm) {

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

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
  // real alpha_J;
  // real theta_J;

}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;
  real<lower=0.1,upper=0.7>  gm;
  real<lower=0.2, upper=0.5> alpha_J;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=340, upper=720> K_CO;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = Ethier_model(Ci[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          alpha_J, theta_J, gm);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // gm ~ normal(0.4, 0.15);
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}
"
stan_Ethier04_D <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real Ethier_model(real Ci,real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J, real gm) {

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

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
  // real alpha_J;
  // real theta_J;

}

parameters {
  real<lower=50, upper=400> V_cmax;
  real<lower=50, upper=400> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;
  real<lower=0.1,upper=0.7>  gm;
  real<lower=0.2, upper=0.5> alpha_J;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=300, upper=900> K_CO;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = Ethier_model(Ci[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          alpha_J, theta_J, gm);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // gm ~ normal(0.4, 0.15);
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}

generated quantities {
  vector[N] log_lik;        // pointwise log-likelihood
  real log_lik_sum;         // total log-likelihood (convenience)

  for (n in 1:N) {
    real A_hat_n;
    A_hat_n = Ethier_model(Ci[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          alpha_J, theta_J, gm);

    // pointwise log-likelihood for Normal(A_hat, sigma)
    log_lik[n] = normal_lpdf(A_obs[n] | A_hat_n, sigma);
  }

  log_lik_sum = sum(log_lik);
}
"
