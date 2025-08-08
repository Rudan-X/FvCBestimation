load_FvCBmodels <- function(){
  models <- list(
    `FvCB (1980)` = list(
      fun  = FvCB1980,
      map  = function(p) list(
        V_cmax=p["Vcmax"], J_max=p["Jmax"], R_d=p["Rd"],
        K_C=p["KC"], K_O=p["KO"], gamma_star=p["Gstar"],
        alpha_J=p["alphaJ"], theta_J=p["thetaJ"], V_tpu=p["Vtpu"]
      )
    ),
    `Harley (1992)` = list(
      fun  = Harley1992,
      map  = function(p) list(
        V_cmax=p["Vcmax"], J_max=p["Jmax"], R_d=p["Rd"],
        K_C=p["KC"], K_O=p["KO"], gamma_star=p["Gstar"], gm=p["gm"]
      )
    ),

    `Ethier & Long (2004)` = list(
      fun  = EthierLong2004,
      map  = function(p) list(
        V_cmax=p["Vcmax"], J_max=p["Jmax"], R_d=p["Rd"],
        K_C=p["KC"], K_O=p["KO"], gamma_star=p["Gstar"],
        gm=p["gm"], alpha_J=p["alphaJ"], theta_J=p["thetaJ"]
      )
    ),
    `von Caemmerer (2000)` = list(
      fun  = Caemmerer2000,
      map  = function(p) list(
        V_cmax=p["Vcmax"], J_max=p["Jmax"], R_d=p["Rd"],
        K_C=p["KC"], K_O=p["KO"], gamma_star=p["Gstar"],
        alpha_J=p["alphaJ"], theta_J=p["thetaJ"],
        V_tpu=p["Vtpu"], alpha_tpu=p["atpu"]
      )
    ),
    `Yin (2004)` = list(
      fun  = Yin2004,
      map  = function(p) {
        phi2m <- if (!is.na(p["phi2"])) p["phi2"] else 0.85
        list(
          V_cmax=p["Vcmax"], J_max=p["Jmax"], R_d=p["Rd"],
          K_C=p["KC"], K_O=p["KO"], gamma_star=p["Gstar"],
          theta_J=p["thetaJ"], phi2m = phi2m,
          f_Q = p["fQ"], f_pseudo = p["fpseudo"], f_cyc = p["fcyc"], h = p["h"]
        )
      }
    ),
    `Dubois (2007)` = list(
      fun  = Dubois2007,
      map  = function(p) list(
        V_cmax=p["Vcmax"], J_max=p["J"], R_d=p["Rd"],
        K_C=p["KC"], K_O=p["KO"], gamma_star=p["Gstar"]
      )
    )
  )
  return(models)
}

