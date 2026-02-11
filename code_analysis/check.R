call_function_FvCB80 <- function(Ci, O, PPFD, pars) {
  model_FvCB80(Ci, 0L, 0, O, PPFD,
                  pars$V_cmax, pars$J_max, pars$R_d,
                  pars$K_C, pars$K_O, pars$gamma_star,
                  pars$alpha_J, pars$theta_J)
}

call_function_FvCB80(c(100,300), 210, 1800, pars)

call_function_Harley92 <- function(Ci, O, PPFD, pars) {
  model_Harley92(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
                    pars$alpha_J,pars$gm)
}
call_function_Harley92(c(100,300), 210, 1800, pars)



pars = list(
  
  # Fixed values:
  g_wp = 1.35,
  g_ch = 0.34,
  
  # estimated values: (focused on Rubisco limited phase only)
  V_cmax=35,
  R_d=0.54,
  K_C=259,
  K_O=179,
  S_co=2592,
  gamma_star=37,
  # Added
  J_max=1.75*35,
  # added to have difference in PPFD response
  alpha_J=0.385, # not specified, use FvCB
  theta_J=0.56 # Ögren and Evans 1993
  
)
call_function_Tholen12 <- function(Ci, O, PPFD, pars) {
  model_Tholen12(Ci, O, PPFD,
                    pars$V_cmax, pars$J_max, pars$R_d,
                    pars$K_C, pars$K_O, pars$gamma_star,
                    pars$alpha_J,pars$theta_J, pars$g_ch ,pars$g_wp)
}
call_function_Tholen12(c(100,300), c(210,210), c(1800,1800), pars)
