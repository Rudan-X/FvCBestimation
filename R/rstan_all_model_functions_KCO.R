stan_all_functions_KCO <- "
functions {
  real J_nrh(real PPFD, real alpha_J, real theta_J, real J_max) {
    real b = -(alpha_J * PPFD + J_max);
    real disc = fmax(0, b*b - 4 * theta_J * alpha_J * PPFD * J_max);
    return (-b - sqrt(disc)) / (2 * theta_J);
  }

  real J_rh(real PPFD, real alpha_J, real J_max) {
    real J  = alpha_J * PPFD / sqrt(1 + square(alpha_J * PPFD / J_max));
    return J;
  }
  real Wc_from(real V_cmax, real C, real Kco) { return V_cmax * C / (C + Kco); }
  real Wj_from(real J, real C, real gamma_star) { return J * C / (4*C + 8*gamma_star); }

  real Wp_from(real V_tpu, real C, real gamma_star, real alpha_G) {
    real denom = C - (1 + 3*alpha_G) * gamma_star;
    return (denom <= 0) ? positive_infinity() : 3 * V_tpu * C / denom;
  }
  real Cref(real Cc, real Ci, int has_Cc) { return has_Cc==1 ? Cc : Ci; }

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

  real function_FvCB80(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J) {

    real Cstar = Cref(Cc, Ci, has_Cc);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Wc    = Wc_from(V_cmax, Cstar, K_CO);
    real Wj    = Wj_from(J, Cstar, gamma_star);
    real Wmin  = fmin(Wc, Wj);
    return (1 - gamma_star / Cstar) * Wmin - R_d;
    }

    real function_Harley92(real Ci, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real gm) {

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


    real function_Tholen12(real Ci, real O, real PPFD,
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

    real function_Caemmerer00(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J, real V_tpu, real alpha_G) {

    real Cstar = Cref(Cc, Ci, has_Cc);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Wc    = Wc_from(V_cmax, Cstar, K_CO);
    real Wj    = Wj_from(J, Cstar, gamma_star);
    real Wp    = Wp_from(V_tpu, Cstar, gamma_star, alpha_G);
    real Wmin  = fmin(Wc, fmin(Wj, Wp));

    return (1 - gamma_star / Cstar) * Wmin - R_d;
    }

    real function_Ethier04(real Ci, real O, real PPFD,
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


    real function_Dubois07(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real alpha_J, real theta_J, real V_tpu, real alpha_G) {

    real Cstar = Cref(Cc, Ci, has_Cc);

    real J = J_nrh(PPFD, alpha_J, theta_J, J_max);

    real Ac    = (1 - gamma_star / Cstar) * Wc_from(V_cmax, Cstar, K_CO) - R_d;
    real Aj    = (1 - gamma_star / Cstar) * Wj_from(J, Cstar, gamma_star) - R_d;
    real Ap    = (1 - gamma_star / Cstar) * Wp_from(V_tpu, Cstar, gamma_star, alpha_G) - R_d;

    return fmin(Ac, fmin(Aj, Ap));
    }



    real function_Yin04(real Ci,int has_Cc, real Cc, real O, real PPFD,
                  real V_cmax, real J_max, real R_d,
                  real K_CO, real gamma_star,
                  real phi2m, real theta_J, real f_Q, real f_pseudo, real h) {

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

  real function_Busch18(real Ci, int has_Cc, real Cc,
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

    vector function_Xiao20(real Ci, real O, real PPFD,
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

    real An = fmin(Ac, Aj);

    real Cc = Ci - An/gm;
    real J_a = ((An +  R_d)*(4*Cc + 8*gamma_star))/(Cc - gamma_star);
    real J_f = fmin(J_a, J_max);
    real Phi2 = J_f/(PPFD * s);

    vector[2] out;
    out[1] = An;
    out[2] = Phi2;
    return out;
    }
}
"
