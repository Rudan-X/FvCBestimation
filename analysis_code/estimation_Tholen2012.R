# --- Packages
suppressWarnings({
  if (!requireNamespace("minpack.lm", quietly = TRUE)) install.packages("minpack.lm")
})
library(minpack.lm)

# ---------- Constants / "true" parameters (Table 1 + text)
true <- list(
  Vcmax = 35.0,          # mmol m^-2 s^-1
  Kc    = 259.0,         # mbar
  Ko    = 179.0,         # mbar
  Rd    = 0.54,          # mmol m^-2 s^-1
  Sco   = 2592.0         # bar bar^-1  (tau = Sc/o)
)

# Light-saturated (Rubisco-limited) branch only; Γ* from τ via Γ* = O/(2τ)
Gamma_star <- function(O, Sco) O / (2 * Sco)          # mbar

# Rubisco-limited carboxylation rate Wc(Cc)
Wc_fun <- function(Cc, Vcmax, Kc, Ko, O) {
  Vcmax * (Cc) / (Cc + Kc * (1 + O / Ko))
}

# Net A given Cc and parameters (Rubisco-limited, standard identity)
A_from_Cc <- function(Cc, Vcmax, Kc, Ko, O, Rd, Sco) {
  Wc <- Wc_fun(Cc, Vcmax, Kc, Ko, O)
  gamma_s <- Gamma_star(O, Sco)
  Wc * (1 - gamma_s / Cc) - Rd
}

# Photorespiratory flux F (needed for variable-gm closure)
F_from_Cc <- function(Cc, Vcmax, Kc, Ko, O, Sco) {
  Wc <- Wc_fun(Cc, Vcmax, Kc, Ko, O)
  gamma_s <- Gamma_star(O, Sco)
  (gamma_s * Wc) / Cc
}

# ---------- Closures to solve Cc

# (1) Constant-gm closure: Ci - Cc = A/gm
solve_Cc_const_gm <- function(Ci, O, pars, gm) {
  fn <- function(Cc) Ci - Cc - A_from_Cc(Cc, pars$Vcmax, pars$Kc, pars$Ko, O, pars$Rd, pars$Sco)/gm
  lo <- 1e-6; hi <- max(Ci, 1.2 * Ci + 5)
  if (fn(lo) * fn(hi) > 0) hi <- hi + 100
  uniroot(fn, c(lo, hi), tol = 1e-9)$root
}

# (2) Variable-gm (Tholen) closure: Ci - Cc = A/g_wp + (F+Rd)/g_ch   (Eq. 6 & 9)
solve_Cc_variable_gm <- function(Ci, O, pars, g_wp, g_ch) {
  fn <- function(Cc) {
    A  <- A_from_Cc(Cc, pars$Vcmax, pars$Kc, pars$Ko, O, pars$Rd, pars$Sco)
    F  <- F_from_Cc(Cc, pars$Vcmax, pars$Kc, pars$Ko, O, pars$Sco)
    Ci - Cc - (A / g_wp + (F + pars$Rd) / g_ch)
  }
  lo <- 1e-6; hi <- max(Ci, 1.2 * Ci + 5)
  # widen if needed
  f_lo <- fn(lo); f_hi <- fn(hi)
  if (f_lo * f_hi > 0) hi <- hi + 200
  uniroot(fn, c(lo, hi), tol = 1e-9)$root
}

# ---------- Simulate datasets

simulate_dataset <- function(Ci_vec, O_vec, pars, mode = c("constant_gm","variable_gm"),
                             gm_const = 0.27, g_wp = 1.35, g_ch = 0.34) {
  mode <- match.arg(mode)
  out <- data.frame()
  for (O in O_vec) {
    for (Ci in Ci_vec) {
      Cc <- switch(mode,
                   constant_gm = solve_Cc_const_gm(Ci, O, pars, gm_const),
                   variable_gm = solve_Cc_variable_gm(Ci, O, pars, g_wp, g_ch)
      )
      A  <- A_from_Cc(Cc, pars$Vcmax, pars$Kc, pars$Ko, O, pars$Rd, pars$Sco)
      out <- rbind(out, data.frame(Ci = Ci, O = O, Cc = Cc, A = A))
    }
  }
  out$mode <- mode
  out
}

# p_i grid (replace with the exact von Caemmerer 1994 grid if you have it)
Ci_grid <- c(20, 30, 40, 50, 60, 80, 100, 120, 150, 200, 250)  # mbar
# Oxygen levels used to recover Kc, Ko, Sco: low and ambient
O_levels <- c( 210, 10)  # mbar (≈1% and 21% O2)

set.seed(1)
dat_const <- simulate_dataset(Ci_grid, O_levels, true, "constant_gm", gm_const = 0.27)
dat_var   <- simulate_dataset(Ci_grid, O_levels, true, "variable_gm", g_wp = 1.35, g_ch = 0.34)

# ---------- Fitting assuming constant gm (Ethier & Livingston style)

# Model wrapper: predict A for given params under constant-gm closure
predict_A_const_gm <- function(par, df, gm_fix = 0.27) {
  Vcmax <- par["Vcmax"]; Kc <- par["Kc"]; Ko <- par["Ko"]; Rd <- par["Rd"]; Sco <- par["Sco"]
  pars <- list(Vcmax = Vcmax, Kc = Kc, Ko = Ko, Rd = Rd, Sco = Sco)
  vapply(seq_len(nrow(df)), function(i) {
    Cc <- solve_Cc_const_gm(df$Ci[i], df$O[i], pars, gm_fix)
    A_from_Cc(Cc, Vcmax, Kc, Ko, df$O[i], Rd, Sco)
  }, numeric(1))
}

fit_constant_gm <- function(df, gm_fix = 0.27, start = NULL, lower = NULL, upper = NULL) {
  if (is.null(start)) {
    start <- c(Vcmax = 20, Kc = 100, Ko = 100, Rd = 0.1, Sco = 2000)
  }
  if (is.null(lower)) {
    lower <- c(Vcmax =  1, Kc =  50, Ko =  50, Rd = 0.0, Sco = 1500)
  }
  if (is.null(upper)) {
    upper <- c(Vcmax = 80, Kc = 500, Ko = 400, Rd = 5.0, Sco = 5000)
  }
  fn <- function(par) predict_A_const_gm(par, df, gm_fix = gm_fix) - df$A
  res <- minpack.lm::nls.lm(par = start, fn = fn, lower = lower, upper = upper,
                            control = minpack.lm::nls.lm.control(maxiter = 300))
  est <- coef(res)
  list(est = est, conv = res$info, message = res$message)
}

fit_const  <- fit_constant_gm(dat_const, gm_fix = 0.27)
fit_var    <- fit_constant_gm(dat_var,   gm_fix = 0.27)

round(fit_const$est, 3)
round(fit_var$est, 3)

# Tidy a compact “Table 1” view
tab1 <- data.frame(
  Parameter = c("Vcmax (mmol m^-2 s^-1)", "Kc (mbar)", "Ko (mbar)", "Rd (mmol m^-2 s^-1)", "Sc/o (bar bar^-1)"),
  Simulated_constant_gm = c(fit_const$est["Vcmax"], fit_const$est["Kc"], fit_const$est["Ko"], fit_const$est["Rd"], fit_const$est["Sco"]),
  Simulated_variable_gm = c(fit_var$est["Vcmax"],   fit_var$est["Kc"],   fit_var$est["Ko"],   fit_var$est["Rd"],   fit_var$est["Sco"])
)
tab1$`Change_vs_constant(%)` <- 100 * (tab1$Simulated_variable_gm - tab1$Simulated_constant_gm) / tab1$Simulated_constant_gm


print(round(tab1, 2))
