library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)

devtools::load_all()

Ci     <- c(50,100,150,200,300,500,600,700,800,1000,1200)
PPFDs  <- c(50,150,300,500,800,1100,1500,1800)

envs_template <- list(
  C_i = Ci,
  O = rep(210, length(Ci))
)


pars <- list(V_cmax = 104,
          J_max = 200,
          R_d = 1.29,
          K_C = 270,
          K_O = 165,
          gamma_star = 34.2,
          alpha_J = 0.3,
          theta_J = 0.8,
          V_tpu = 11.5,
          gm = 0.2,
          alpha_G = 0.3,
          phi2m = 0.85,
          f_Q = 1,
          f_pseudo = 0,
          f_cyc = 1000,
          h = 14/3,
          g_wp = 2.35,
          g_ch = 0.34,
          max_alpha_G = 0.09,
          max_alpha_S = 0.38,
          N_max = 1.21,
          s=0.369,
          Phi2LL = 0.75
         )

FvCB98 <- model_FvCB1980()
Caemmerer00 <- model_Caemmerer2000(base = FvCB98, ci_tpu_min = 100)
Ethier04 <- model_Ethier2004()
models <- list(FvCB98,model_Harley1992(base = FvCB98),
               Caemmerer00,Ethier04,
               model_Yin2004(base = FvCB98),
               model_Dubois2007(base = Caemmerer00, ci_tpu_min = 100),
               model_Tholen2012(base = FvCB98),model_Busch2018(base = FvCB98, ci_tpu_min = 100),
               model_Xiao2021(base = Ethier04))


model_names <- c("FvCB (1980) (Wc,Wj) (Ci)", "Harley (1992) (Ac,Aj) (Cc)",
                 "von Caemmerer (2000) (Wc,Wj,Wp) (Ci)","Ethier & Livingston (2004) (Ac,Aj) (Cc)",
                 "Yin (2004) (Wc,Wj) (Ci)","Dubois (2007) (Ac,Aj,Ap) (Ci)",
                 "Tholen (2012) (Wc,Wj) (Cc)","Busch (2018) (Wc,Wj,Wp) (Ci)",
                 "Xiao (2021) (Ac,Aj) (Cc)")

names(models) <- model_names
run_model_once <- function(model_obj, envs_base, pars, PPFD_val) {
  n <- length(envs_base$C_i)
  envs <- list(
    C_i  = envs_base$C_i,
    O    = rep(envs_base$O, n),     # ensure same length as C_i
    PPFD = rep(PPFD_val, n)
  )
  res  <- predict_A(model_obj, envs, pars)

  # if (model_obj$name=="Xiao2021"){
  #   temp <- predict_Phi2(model_obj,res$An,envs,pars)
  #   res <- cbind(res, temp)
  # }
  tibble(Ci = envs$C_i, A = res$An, PPFD = PPFD_val)
}

df <- imap_dfr(models, function(model_obj, model_name) {
  map_dfr(PPFDs, ~ run_model_once(model_obj, envs_template, pars, .x)) %>%
    mutate(Model = model_name)
})


df$PPFD <- paste0("PPFD: ", df$PPFD)
df$PPFD <-   factor(df$PPFD, levels=unique(df$PPFD))
# df$Model <- factor(df$Model, levels= new_names)

ggplot(df, aes(Ci, A, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  facet_wrap(~PPFD, ncol = 4, dir = "v") +
  labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
       y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
       color = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("results/Figures/Figure1.b_Single_parameter.png",width = 12, height = 7)
