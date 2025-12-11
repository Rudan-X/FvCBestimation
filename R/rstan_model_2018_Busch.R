stan_Busch18 <- "
functions{
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }

  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }


  // ---------- Busch 2018 alpha helpers ----------
  real alpha_beta(real max_aG, real max_aS) {
      return 3*max_aG / (3*max_aG + 2*max_aS);
    }
  real phi_from(real Cstar, real gamma_star) {
    return 2 * gamma_star / fmax(Cstar, 1e-9);
  }
  real alpha_G_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = N_max * beta / fmax(Vo_c, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = (3 * N_max * (1 - beta)) / (2 * fmax(Vo_c, 1e-12));
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aG;
    real val = 4 * N_max * beta * (1/phi + 1) / (J - thresh);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aS;
    real val = 6 * N_max * (1 - beta) * (1/phi + 1) / (J - thresh);
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = N_max * beta * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = 1.5 * N_max * (1 - beta) * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aS, fmax(0, val));
  }

  real model_Busch2018(real Ci, int has_Cc, real Cc,
                    real O, real PPFD,
                    real Vcmax, real J_max, real R_d,
                    real K_C, real K_O, real gamma_star,
                    real alpha_J, real theta_J, real Vtpu,
                    real N_max, real max_aG, real max_aS) {

    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real K_CO = K_C * (1 + O_conv / Ko_conv);
    real Cstar = Cref(Cc, Ci, has_Cc);

    real J     = J_nrh(PPFD, alpha_J, theta_J, J_max);
    real Wc    = Wc_from(Vcmax, Cstar, K_CO);
    real aGj   = alpha_G_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real aSj   = alpha_S_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real phi   = phi_from(Cstar, gamma_star);
    real Wj    = J / (4 + (4 + 8*aGj + 4*aSj) * phi);
    real aGp   = alpha_G_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real aSp   = alpha_S_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real denom = 1 - 0.5 * (1 + 3*aGp + 4*aSp) * phi;
    real Wp    = (denom <= 0) ? positive_infinity() : 3 * Vtpu / denom;
    // choose limitation (piecewise is okay; if HMC struggles, switch to softmin)
    real Wmin  = fmin(Wc, fmin(Wj, Wp));
    // pick matching gamma* for the limiting regime:

    real g_c = (1 - alpha_G_c(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS)) * gamma_star;
    real g_j = (1 - aGj) * gamma_star;
    real g_p = (1 - aGp) * gamma_star;
    real g_pick = (Wmin == Wc) ? g_c : ((Wmin == Wj) ? g_j : g_p);
    return (1 - g_pick / Cstar) * Wmin - R_d;
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
  real<lower=8, upper=20> V_tpu;
  real<lower=0.2, upper=0.5> alpha_J;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=200, upper=300> K_C;
  real<lower=150, upper=300> K_O;
  real<lower=0.01, upper=0.2> max_aG;
  real<lower=0.01, upper=0.5> max_aS;
  real<lower=0.1, upper=2> N_max;
  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = model_Busch2018(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_C, K_O, gamma_star,
                          alpha_J, theta_J, V_tpu,
                          N_max, max_aG, max_aS);
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


stan_Busch18_B <- "
functions{
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }

  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }


  // ---------- Busch 2018 alpha helpers ----------
  real alpha_beta(real max_aG, real max_aS) {
      return 3*max_aG / (3*max_aG + 2*max_aS);
    }
  real phi_from(real Cstar, real gamma_star) {
    return 2 * gamma_star / fmax(Cstar, 1e-9);
  }
  real alpha_G_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = N_max * beta / fmax(Vo_c, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = (3 * N_max * (1 - beta)) / (2 * fmax(Vo_c, 1e-12));
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aG;
    real val = 4 * N_max * beta * (1/phi + 1) / (J - thresh);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aS;
    real val = 6 * N_max * (1 - beta) * (1/phi + 1) / (J - thresh);
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = N_max * beta * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = 1.5 * N_max * (1 - beta) * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aS, fmax(0, val));
  }

  real model_Busch2018(real Ci, int has_Cc, real Cc,
                    real O, real PPFD,
                    real Vcmax, real J_max, real R_d,
                    real K_C, real K_O, real gamma_star,
                    real alpha_J, real theta_J, real Vtpu,
                    real N_max, real max_aG, real max_aS) {

    real O_conv  = O * 1000;
    real Ko_conv = K_O * 1000;

    real K_CO = K_C * (1 + O_conv / Ko_conv);
    real Cstar = Cref(Cc, Ci, has_Cc);

    real J     = J_nrh(PPFD, alpha_J, theta_J, J_max);
    real Wc    = Wc_from(Vcmax, Cstar, K_CO);
    real aGj   = alpha_G_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real aSj   = alpha_S_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real phi   = phi_from(Cstar, gamma_star);
    real Wj    = J / (4 + (4 + 8*aGj + 4*aSj) * phi);
    real aGp   = alpha_G_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real aSp   = alpha_S_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);

    real denom = 1 - 0.5 * (1 + 3*aGp + 4*aSp) * phi;
    real Wp;
    if (denom <= 0) {
      Wp = 1e6;                 // large positive -> will not limit
    } else {
      Wp = 3 * Vtpu / denom;
    }

    // choose limitation (piecewise is okay; if HMC struggles, switch to softmin)
    real Wmin  = fmin(Wc, fmin(Wj, Wp));
    // pick matching gamma* for the limiting regime:

    real g_c = (1 - alpha_G_c(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS)) * gamma_star;
    real g_j = (1 - aGj) * gamma_star;
    real g_p = (1 - aGp) * gamma_star;
    real g_pick = (Wmin == Wc) ? g_c : ((Wmin == Wj) ? g_j : g_p);
    return (1 - g_pick / Cstar) * Wmin - R_d;
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
  real alpha_G;
  real V_tpu;
  real max_aG;
  real max_aS;
  real N_max;

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
    A_hat[n] = model_Busch2018(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_C, K_O, gamma_star,
                          alpha_J, theta_J, V_tpu,
                          N_max, max_aG, max_aS);
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


stan_Busch18_C <- "
functions{
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }

  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }


  // ---------- Busch 2018 alpha helpers ----------
  real alpha_beta(real max_aG, real max_aS) {
      return 3*max_aG / (3*max_aG + 2*max_aS);
    }
  real phi_from(real Cstar, real gamma_star) {
    return 2 * gamma_star / fmax(Cstar, 1e-9);
  }
  real alpha_G_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = N_max * beta / fmax(Vo_c, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = (3 * N_max * (1 - beta)) / (2 * fmax(Vo_c, 1e-12));
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aG;
    real val = 4 * N_max * beta * (1/phi + 1) / (J - thresh);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aS;
    real val = 6 * N_max * (1 - beta) * (1/phi + 1) / (J - thresh);
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = N_max * beta * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = 1.5 * N_max * (1 - beta) * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aS, fmax(0, val));
  }

  real model_Busch2018(real Ci, int has_Cc, real Cc,
                    real O, real PPFD,
                    real Vcmax, real J_max, real R_d,
                    real K_CO, real gamma_star,
                    real alpha_J, real theta_J, real Vtpu,
                    real N_max, real max_aG, real max_aS) {

    real Cstar = Cref(Cc, Ci, has_Cc);

    real J     = J_nrh(PPFD, alpha_J, theta_J, J_max);
    real Wc    = Wc_from(Vcmax, Cstar, K_CO);
    real aGj   = alpha_G_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real aSj   = alpha_S_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real phi   = phi_from(Cstar, gamma_star);
    real Wj    = J / (4 + (4 + 8*aGj + 4*aSj) * phi);
    real aGp   = alpha_G_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real aSp   = alpha_S_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real denom = 1 - 0.5 * (1 + 3*aGp + 4*aSp) * phi;
    real Wp    = (denom <= 0) ? positive_infinity() : 3 * Vtpu / denom;
    // choose limitation (piecewise is okay; if HMC struggles, switch to softmin)
    real Wmin  = fmin(Wc, fmin(Wj, Wp));
    // pick matching gamma* for the limiting regime:

    real g_c = (1 - alpha_G_c(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS)) * gamma_star;
    real g_j = (1 - aGj) * gamma_star;
    real g_p = (1 - aGp) * gamma_star;
    real g_pick = (Wmin == Wc) ? g_c : ((Wmin == Wj) ? g_j : g_p);
    return (1 - g_pick / Cstar) * Wmin - R_d;
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
  real<lower=8, upper=20> V_tpu;
  real<lower=0.2, upper=0.5> alpha_J;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=340, upper=720> K_CO;
  real<lower=0.01, upper=0.2> max_aG;
  real<lower=0.01, upper=0.5> max_aS;
  real<lower=0.1, upper=2> N_max;
  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = model_Busch2018(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          alpha_J, theta_J, V_tpu,
                          N_max, max_aG, max_aS);
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


stan_Busch18_D <- "
functions{
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }
  real Wc_from(real V_cmax, real C, real K_CO) { return V_cmax * C / (C + K_CO); }

  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }


  // ---------- Busch 2018 alpha helpers ----------
  real alpha_beta(real max_aG, real max_aS) {
      return 3*max_aG / (3*max_aG + 2*max_aS);
    }
  real phi_from(real Cstar, real gamma_star) {
    return 2 * gamma_star / fmax(Cstar, 1e-9);
  }
  real alpha_G_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = N_max * beta / fmax(Vo_c, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_c(real Cstar, real Vcmax, real K_CO, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real Vo_c = Vcmax * (2*gamma_star) / (Cstar + K_CO);
    real beta = alpha_beta(max_aG, max_aS);
    real val  = (3 * N_max * (1 - beta)) / (2 * fmax(Vo_c, 1e-12));
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aG;
    real val = 4 * N_max * beta * (1/phi + 1) / (J - thresh);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_j(real Cstar, real J, real gamma_star,
                 real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real thresh = N_max * (2*beta + 6);
    if (J <= thresh) return max_aS;
    real val = 6 * N_max * (1 - beta) * (1/phi + 1) / (J - thresh);
    return fmin(max_aS, fmax(0, val));
  }
  real alpha_G_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = N_max * beta * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aG, fmax(0, val));
  }
  real alpha_S_p(real Cstar, real gamma_star,
                 real Vtpu, real N_max, real max_aG, real max_aS) {
    real beta = alpha_beta(max_aG, max_aS);
    real phi  = phi_from(Cstar, gamma_star);
    real num  = 1.5 * N_max * (1 - beta) * fmax(2/phi - 1, 0);
    real den  = 6 * Vtpu + 3 * N_max * (2 - beta);
    real val  = num / fmax(den, 1e-12);
    return fmin(max_aS, fmax(0, val));
  }

  real model_Busch2018(real Ci, int has_Cc, real Cc,
                    real O, real PPFD,
                    real Vcmax, real J_max, real R_d,
                    real K_CO, real gamma_star,
                    real alpha_J, real theta_J, real Vtpu,
                    real N_max, real max_aG, real max_aS) {

    real Cstar = Cref(Cc, Ci, has_Cc);

    real J     = J_nrh(PPFD, alpha_J, theta_J, J_max);
    real Wc    = Wc_from(Vcmax, Cstar, K_CO);
    real aGj   = alpha_G_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real aSj   = alpha_S_j(Cstar, J, gamma_star, N_max, max_aG, max_aS);
    real phi   = phi_from(Cstar, gamma_star);
    real Wj    = J / (4 + (4 + 8*aGj + 4*aSj) * phi);
    real aGp   = alpha_G_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real aSp   = alpha_S_p(Cstar, gamma_star, Vtpu, N_max, max_aG, max_aS);
    real denom = 1 - 0.5 * (1 + 3*aGp + 4*aSp) * phi;
    real Wp    = (denom <= 0) ? positive_infinity() : 3 * Vtpu / denom;
    // choose limitation (piecewise is okay; if HMC struggles, switch to softmin)
    real Wmin  = fmin(Wc, fmin(Wj, Wp));
    // pick matching gamma* for the limiting regime:

    real g_c = (1 - alpha_G_c(Cstar, Vcmax, K_CO, gamma_star, N_max, max_aG, max_aS)) * gamma_star;
    real g_j = (1 - aGj) * gamma_star;
    real g_p = (1 - aGp) * gamma_star;
    real g_pick = (Wmin == Wc) ? g_c : ((Wmin == Wj) ? g_j : g_p);
    return (1 - g_pick / Cstar) * Wmin - R_d;
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
  real<lower=8, upper=20> V_tpu;
  real<lower=0.2, upper=0.5> alpha_J;
  real<lower=0.3, upper=0.9> theta_J;
  real<lower=300, upper=900> K_CO;
  real<lower=0.01, upper=0.2> max_aG;
  real<lower=0.01, upper=0.5> max_aS;
  real<lower=0.1, upper=2> N_max;
  real<lower=0> sigma;
}

model {
  vector[N] A_hat;

  // Likelihood
  for (n in 1:N) {
    A_hat[n] = model_Busch2018(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          alpha_J, theta_J, V_tpu,
                          N_max, max_aG, max_aS);
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
    A_hat_n = model_Busch2018(Ci[n], has_Cc, Cc[n], O[n], PPFD[n],
                          V_cmax, J_max, R_d,
                          K_CO, gamma_star,
                          alpha_J, theta_J, V_tpu,
                          N_max, max_aG, max_aS);

    // pointwise log-likelihood for Normal(A_hat, sigma)
    log_lik[n] = normal_lpdf(A_obs[n] | A_hat_n, sigma);
  }

  log_lik_sum = sum(log_lik);
}
"


