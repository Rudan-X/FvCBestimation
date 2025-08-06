calc_Anet <- function(Ci, Vcmax, Jmax, gamma_s, KmCO){ #, Rd, KO
  gi <- 0.2
  Rd <- 0.5
  a <- -1 / gi
  b <- ((Vcmax - Rd) / gi) + Ci + KmCO
  c <- Rd * (Ci + KmCO) - Vcmax * (Ci - gamma_s)

  sqrt_term <- sqrt(b^2 - 4 * a * c)
  Ac <- (-b + sqrt_term) / (2 * a)
  return(Ac)
}


C <- c(50,100,150,200,300,400, 600, 800, 1200)
true.val <- c(56,112, 45, 400)
An <- calc_Anet(C, true.val[1], true.val[2], true.val[3],true.val[4])

ACi_curve <- data.frame(Ci = C, A = An)

plot(ACi_curve$Ci, ACi_curve$A)
# Initial guesses
start_vals <- list(Vcmax = 50, gi = 0.001, Rd=0.5)

# Fit model
fit <- nls(
  A ~ ac_quad(Ci, Vcmax, gi, Rd),
  data = ACi_curve,
  start = start_vals,
  algorithm = "port",
  lower = c(Vcmax = 0, Rd = -3),
  upper = c(Vcmax = 500, Rd = 10),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

coef(fit)
summary(fit)

estA <- predict(fit)
plot(ACi_curve$Ci, estA)
