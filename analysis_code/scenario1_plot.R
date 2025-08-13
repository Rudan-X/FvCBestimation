library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

devtools::load_all()

# Assuming 'scenarios' is already computed from the previous simulation code.
# If not, include the simulation code here to generate 'scenarios'.

# Grid search for environmental input (unchanged):
Ci <- c(20,50,100,150,200,300,400,500,600,700,800,1000,1200)
PPFDs <- c(20,100,300,500,700,1000,1200,1500)

# No need for envs_template or param here, as we're using precomputed 'scenarios'

# No need for manual models list, as models are in scenarios

# Function to compute mean A for one model in one scenario



compute_mean_A <- function(an_array, Ci_vals, PPFD_vals) {
  # an_array is 3D [Ci, PPFD, samples]
  mean_A <- apply(an_array, c(1, 2), mean, na.rm = TRUE)
  # Convert to long format tibble
  tibble(
    Ci = rep(Ci_vals, times = length(PPFD_vals)),
    PPFD = rep(PPFD_vals, each = length(Ci_vals)),
    A = as.vector(mean_A)
  )
}

# Loop over scenarios to create and save plots
for (scenario_name in names(scenarios)) {
  # For each scenario, collect data across models
  df_scenario <- imap_dfr(scenarios[[scenario_name]], function(model_res, model_name) {
    # model_res has $An (3D array)
    df_model <- compute_mean_A(model_res$An, Ci, PPFDs) %>%
      mutate(Model = model_name)
  })
  
  df_scenario$PPFD <- factor(df_scenario$PPFD, levels = PPFDs)
  
  # Create the plot (similar to original, but using means)
  p <- ggplot(df_scenario, aes(Ci, A, color = Model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    facet_wrap(~PPFD, ncol = 4) +
    labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
         y = expression(A~"(µmol m"^{-2}~s^{-1}*")"),
         color = "Model",
         title = paste("Mean A for", scenario_name)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  # Save the plot
  plot_file <- paste0("plot_", scenario_name, ".png")
  ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
}
p