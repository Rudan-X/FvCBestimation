stan_Yin04 <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real Kco) { return V_cmax * C / (C + Kco); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Yin04(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real phi2m, real theta_J, real f_Q, real f_pseudo, real h) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real Cstar = Cref(Cc, Ci, has_Cc);
    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real Wc    = Wc_from(V_cmax, Cstar, K_CO);

    real a = 4*Cc + 8*gamma_star;
    real b = 3*Cc + 7*gamma_star;

    // f_cyc = ( a*(2+f_Q)/(h*b) - (1 - f_pseudo) ) / ( a/(h*b) - 1 )
    real denom = a/(h*b) - 1;
    // guard denom near zero
    denom = (fabs(denom) < 1e-12) ? (denom >= 0 ? 1e-12 : -1e-12) : denom;
    real f_cyc = ( a*(2 + f_Q)/(h*b) - (1 - f_pseudo) ) / denom;


    // Clamp the fractions to [0, 1] to keep them physically meaningful
    f_cyc = fmin(1, fmax(0, f_cyc));

    real one_minus_fc = fmax(0, 1 - f_cyc);
    real alpha_J = one_minus_fc / (1 + (one_minus_fc / fmax(phi2m, 1e-12)));

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Wj    = J * ((2 + f_Q - f_cyc) * Cstar) / (h * (3*Cstar + 7*gamma_star) * fmax(1e-12, (1 - f_cyc)));

    real Wmin  = fmin(Wc, Wj);

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  // real alpha_G;

  // flags and model selection
  int<lower=0,upper=1> has_Cc;     // 1 if Cc provided (e.g., Harley variants)
  vector[N] Cc;                    // used only if has_Cc==1
}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;

  real<lower=0.3, upper=0.9> theta_J;
  real<lower=200, upper=300> K_C;
  real<lower=150, upper=300> K_O;

  real<lower=0.5, upper=0.9> phi2m;
  real<lower=0.5, upper=1> f_Q;
  real<lower=0, upper=0.3> f_pseudo;
  real<lower=3, upper=14/3> h;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = model_Yin04(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_C, K_O, gamma_star,
                          phi2m, theta_J, f_Q, f_pseudo, h);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // V_tpu ~ normal(14, 3);   // within 8–20
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}
"

stan_Yin04_B <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real Kco) { return V_cmax * C / (C + Kco); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Yin04(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_C, real K_O, real gamma_star,
                  real phi2m, real theta_J, real f_Q, real f_pseudo, real h) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)
    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real Cstar = Cref(Cc, Ci, has_Cc);
    real K_CO = K_C * (1 + O_conv / Ko_conv);

    real Wc    = Wc_from(V_cmax, Cstar, K_CO);

    real a = 4*Cc + 8*gamma_star;
    real b = 3*Cc + 7*gamma_star;

    // f_cyc = ( a*(2+f_Q)/(h*b) - (1 - f_pseudo) ) / ( a/(h*b) - 1 )
    real denom = a/(h*b) - 1;
    // guard denom near zero
    denom = (fabs(denom) < 1e-12) ? (denom >= 0 ? 1e-12 : -1e-12) : denom;
    real f_cyc = ( a*(2 + f_Q)/(h*b) - (1 - f_pseudo) ) / denom;


    // Clamp the fractions to [0, 1] to keep them physically meaningful
    f_cyc = fmin(1, fmax(0, f_cyc));

    real one_minus_fc = fmax(0, 1 - f_cyc);
    real alpha_J = one_minus_fc / (1 + (one_minus_fc / fmax(phi2m, 1e-12)));

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Wj    = J * ((2 + f_Q - f_cyc) * Cstar) / (h * (3*Cstar + 7*gamma_star) * fmax(1e-12, (1 - f_cyc)));

    real Wmin  = fmin(Wc, Wj);

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  real phi2m;
  real f_Q;
  real f_pseudo;
  real h;

  // flags and model selection
  int<lower=0,upper=1> has_Cc;     // 1 if Cc provided (e.g., Harley variants)
  vector[N] Cc;                    // used only if has_Cc==1
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
    A_hat[n] = model_Yin04(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_C, K_O, gamma_star,
                          phi2m, theta_J, f_Q, f_pseudo, h);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // V_tpu ~ normal(14, 3);   // within 8–20
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}
"

stan_Yin04_C <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real Kco) { return V_cmax * C / (C + Kco); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Yin04(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real phi2m, real theta_J, real f_Q, real f_pseudo, real h) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)

    real Cstar = Cref(Cc, Ci, has_Cc);

    real Wc    = Wc_from(V_cmax, Cstar, K_CO);

    real a = 4*Cc + 8*gamma_star;
    real b = 3*Cc + 7*gamma_star;

    // f_cyc = ( a*(2+f_Q)/(h*b) - (1 - f_pseudo) ) / ( a/(h*b) - 1 )
    real denom = a/(h*b) - 1;
    // guard denom near zero
    denom = (fabs(denom) < 1e-12) ? (denom >= 0 ? 1e-12 : -1e-12) : denom;
    real f_cyc = ( a*(2 + f_Q)/(h*b) - (1 - f_pseudo) ) / denom;


    // Clamp the fractions to [0, 1] to keep them physically meaningful
    f_cyc = fmin(1, fmax(0, f_cyc));

    real one_minus_fc = fmax(0, 1 - f_cyc);
    real alpha_J = one_minus_fc / (1 + (one_minus_fc / fmax(phi2m, 1e-12)));

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Wj    = J * ((2 + f_Q - f_cyc) * Cstar) / (h * (3*Cstar + 7*gamma_star) * fmax(1e-12, (1 - f_cyc)));

    real Wmin  = fmin(Wc, Wj);

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  // real alpha_G;

  // flags and model selection
  int<lower=0,upper=1> has_Cc;     // 1 if Cc provided (e.g., Harley variants)
  vector[N] Cc;                    // used only if has_Cc==1
}

parameters {
  real<lower=50, upper=200> V_cmax;
  real<lower=100, upper=250> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;

  real<lower=0.3, upper=0.9> theta_J;
  real<lower=340, upper=720> K_CO;

  real<lower=0.5, upper=0.9> phi2m;
  real<lower=0.5, upper=1> f_Q;
  real<lower=0, upper=0.3> f_pseudo;
  real<lower=3, upper=14/3> h;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = model_Yin04(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          phi2m, theta_J, f_Q, f_pseudo, h);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // V_tpu ~ normal(14, 3);   // within 8–20
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}
"

stan_Yin04_D <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real Kco) { return V_cmax * C / (C + Kco); }


  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

  real model_Yin04(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real phi2m, real theta_J, real f_Q, real f_pseudo, real h) {
    // Convert O and K_O to micromol mol^-1 (same scaling as your R code)

    real Cstar = Cref(Cc, Ci, has_Cc);

    real Wc    = Wc_from(V_cmax, Cstar, K_CO);

    real a = 4*Cc + 8*gamma_star;
    real b = 3*Cc + 7*gamma_star;

    // f_cyc = ( a*(2+f_Q)/(h*b) - (1 - f_pseudo) ) / ( a/(h*b) - 1 )
    real denom = a/(h*b) - 1;
    // guard denom near zero
    denom = (fabs(denom) < 1e-12) ? (denom >= 0 ? 1e-12 : -1e-12) : denom;
    real f_cyc = ( a*(2 + f_Q)/(h*b) - (1 - f_pseudo) ) / denom;


    // Clamp the fractions to [0, 1] to keep them physically meaningful
    f_cyc = fmin(1, fmax(0, f_cyc));

    real one_minus_fc = fmax(0, 1 - f_cyc);
    real alpha_J = one_minus_fc / (1 + (one_minus_fc / fmax(phi2m, 1e-12)));

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Wj    = J * ((2 + f_Q - f_cyc) * Cstar) / (h * (3*Cstar + 7*gamma_star) * fmax(1e-12, (1 - f_cyc)));

    real Wmin  = fmin(Wc, Wj);

    return (1 - gamma_star / Cstar) * Wmin - R_d;
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
  // real alpha_G;

  // flags and model selection
  int<lower=0,upper=1> has_Cc;     // 1 if Cc provided (e.g., Harley variants)
  vector[N] Cc;                    // used only if has_Cc==1
}

parameters {
  real<lower=50, upper=400> V_cmax;
  real<lower=50, upper=400> J_max;
  real<lower=0.5, upper=5> R_d;
  real<lower=10, upper=60> gamma_star;

  real<lower=0.3, upper=0.9> theta_J;
  real<lower=300, upper=900> K_CO;

  real<lower=0.5, upper=0.9> phi2m;
  real<lower=0.5, upper=1> f_Q;
  real<lower=0, upper=0.3> f_pseudo;
  real<lower=3, upper=14/3> h;

  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = model_Yin04(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          phi2m, theta_J, f_Q, f_pseudo, h);
  }
  A_obs ~ normal(A_hat, sigma);

  // Priors: weakly informative, centered in the middle of the bounds
  // V_cmax ~ normal(125, 30);       // within 50–200
  // J_max ~ normal(175, 40);        // within 100–250
  // R_d ~ normal(2.5, 1.0);         // within 0.5–5
  // gamma_star ~ normal(35, 10);   // within 10–60
  // V_tpu ~ normal(14, 3);   // within 8–20
  sigma ~ normal(0, 2) T[0,];    // positive residual SD
}

generated quantities {
  vector[N] log_lik;        // pointwise log-likelihood
  real log_lik_sum;         // total log-likelihood (convenience)

  for (n in 1:N) {
    real A_hat_n;
    A_hat_n = model_Yin04(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          phi2m, theta_J, f_Q, f_pseudo, h);

    // pointwise log-likelihood for Normal(A_hat, sigma)
    log_lik[n] = normal_lpdf(A_obs[n] | A_hat_n, sigma);
  }

  log_lik_sum = sum(log_lik);
}
"
