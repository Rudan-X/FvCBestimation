library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)



filen<-paste0("R/")
myls <- list.files(path=filen,pattern="*.R")
myls<-paste0(filen,myls)
invisible(sapply(myls,FUN=source))



devtools::load_all()

# load(file="results/Rdata/fig2_sampled_params.RData")

load(file="results/Rdata/fig2_simulated_curves.RData")

rst_FvCB <- stan_model(model_code = stan_FvCB80)
rst_Harley <- stan_model(model_code = stan_Harley92)
rst_Caemmerer <- stan_model(model_code = stan_Caemmerer00)
rst_Ethier <- stan_model(model_code = stan_Ethier04)
rst_Yin <- stan_model(model_code = stan_Yin04)
rst_Dubois <- stan_model(model_code = stan_Dubois07)
rst_Tholen <- stan_model(model_code = stan_Tholen12)
rst_Busch <- stan_model(model_code = stan_Busch18)
rst_Xiao_A <- stan_model(model_code = stan_Xiao20_An)
rst_Xiao_Phi <- stan_model(model_code = stan_Xiao20_Phi)

rst_models <- list(
  FvCB80 = rst_FvCB,
  Harley92 = rst_Harley,
  Caemmerer00  = rst_Caemmerer,
  Ethier04  = rst_Ethier,
  Yin04 = rst_Yin,
  Dubois07 = rst_Dubois,
  Tholen12 = rst_Tholen,
  Busch18 = rst_Busch,
  Xiao20 = list(rst_Xiao_A,rst_Xiao_Phi)
)


param_candidates <- c("V_cmax","J_max","R_d","gamma_star",
                      "theta_J","alpha_J", "K_C", "K_O","V_tpu", "alpha_G","gm",
                      "phi2m", "f_Q", "f_pseudo", "h", "N_max", "max_aG", "max_aS",
                      "g_ch","g_wp","s","Phi2LL")


build_stan_data <- function(sim_res) {
  sim_res <- sim_res[complete.cases(sim_res),]
  n=length(sim_res$Ci)
  data <- list(
    N     = as.integer(n),
    Ci    = as.numeric(sim_res$Ci),
    O     = as.numeric(sim_res$O),
    PPFD  = as.numeric(sim_res$PPFD),
    A_obs = as.numeric(sim_res$An),
    has_Cc     = 0L,
    Cc         = rep(0.0, n)
  )
  if ("Phi2"%in%colnames(sim_res)){
    data$Phi2_obs <- as.numeric(sim_res$Phi2)
  }
  return(data)
}



datatypes <- c("ACi(satQ)", "ACi & AQ", "ACi(different Q)", "ACi(different Q) & CF")


fit_one_repitition <- function(rep_id,rst_models,   # list of compiled Stan models
                               iter = 3000, warmup = 1000, chains = 4,
                               adapt_delta = 0.95, max_treedepth = 12) {

  res_list <- list()
  # storage
  for (m in 1:length(rst_models)) {

    if (m==9){
      dl <- 4
    }else{dl <- 3}

    for (d in 1:dl){
      if (m==9){
        if (d!=4){
          stan_m <- rst_models[[m]][[1]]
        }else{
          stan_m <- rst_models[[m]][[2]]
        }
      }else{
        stan_m <- rst_models[[m]]
      }

      print(paste0("rep",rep_x,"_",names(rst_models)[m],"_",datatypes[d]))

      case_name <- paste0("rep",rep_id,"_",names(rst_models)[m],"_",datatypes[d])
      x <- which(names(simulated_curves)==case_name)
      sim_res <- simulated_curves[[x]]
      sim_res$O <- 210

      sdata <- build_stan_data(sim_res)
      fit <- sampling(
        stan_m, data = sdata,
        chains = chains, cores=4,iter = iter, warmup = warmup,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        refresh = 0
      )

      dr <- posterior::as_draws_df(fit)
      keep <- intersect(colnames(dr), param_candidates)
      keep <- c(keep, "lp__",".draw")
      dr <- dr[, keep, drop = FALSE]

      res_list[[length(res_list) + 1]] <- dr
      names(res_list)[length(res_list)] <- paste0("rep",rep_x,"_",names(rst_models)[m],"_",datatypes[d])
    }
  }
  save(res_list, file=paste0("results/Rdata/fig1_identifiability_i",rep_id,".RData"))
}

for (rep_x in 8:100){
  fit_one_repitition(rep_x, rst_models = rst_models)
}



res_list$scenario <- factor(res_list$scenario, levels = datatypes[1:3])

common_param <- c("V_cmax","J_max","R_d","gamma_star",
 "K_C", "K_O", "theta_J","alpha_J")


df_temp <- res_list[res_list$param%in%common_param,]
df_temp$param <- factor(df_temp$param, levels = common_param)
df_temp$model <- factor(df_temp$model, levels = names(rst_models))
df_temp$scenario <- factor(df_temp$scenario, levels = c("ACi(satQ)", "ACi & AQ", "ACi(different Q)"))


df_temp2 <- true_params[true_params$param%in%common_param,]

df_temp2$param <- factor(df_temp2$param, levels = common_param)

ggplot(data = df_temp,
       aes(x = scenario, y = value, fill = scenario, colour = scenario)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.30) +
  # dashed truth line per facet (uses df_temp2$truth)
  geom_hline(data = df_temp2,
             aes(yintercept = value),
             inherit.aes = FALSE,
             linetype = "dashed", linewidth = 0.45, alpha = 0.9) +
  facet_grid(param ~ model, scales = "free_y") +
  labs(x = "Scenario", y = NULL, colour = "Scenario", fill = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.spacing.x = grid::unit(8, "pt"),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank()
  )


ggplot() +
  geom_density(data= df_temp,
               aes(x = value,
                   colour   = scenario# ,                        # same hue per replicate
               ), #, linetype = scenario
               # group    = interaction(rep, scenario)
               size = 0.9, adjust = 1) +
  geom_vline(data = df_temp2,
             aes(xintercept = value),
             linetype = "dashed", linewidth = 0.45, alpha = 0.9) +
  facet_wrap(param ~ model, scales = "free", ncol(rst_models)) +
  scale_linetype_manual(values = c(single = "solid", multi = "dashed")) +
  labs(x = NULL, y = "Density",
       colour = "Replicate", linetype = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.spacing.x = grid::unit(8, "pt"))

ggsave(filename = paste0("results/Figures/Fig2_identifiability_models_vs_datatype_i",rep_x,".png"),width = 7, height = 10)

############## Non essential #############
diff_param <- setdiff(param_candidates,common_param)
df_temp <- res_list[res_list$param%in%diff_param,]
df_temp$param <- factor(df_temp$param, levels = diff_param)
df_temp$model <- factor(df_temp$model, levels = names(rst_models))
df_temp$scenario <- factor(df_temp$scenario, levels = c("ACi(satQ)", "ACi & AQ", "ACi(different Q)"))


df_temp2 <- true_params[true_params$param%in%diff_param,]

df_temp2$param <- factor(df_temp2$param, levels = diff_param)


ggplot(data = df_temp,
       aes(x = scenario, y = value, fill = scenario, colour = scenario)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.30) +
  # dashed truth line per facet (uses df_temp2$truth)
  geom_hline(data = df_temp2,
             aes(yintercept = value),
             inherit.aes = FALSE,
             linetype = "dashed", linewidth = 0.45, alpha = 0.9) +
  facet_grid(param ~ model, scales = "free_y") +
  labs(x = "Scenario", y = NULL, colour = "Scenario", fill = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.spacing.x = grid::unit(8, "pt"),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank()
  )

ggsave(filename = paste0("results/Figures/Fig2_identifiability_models_vs_datatype_specific_i",rep_x,".png"),width = 7, height = 10)

