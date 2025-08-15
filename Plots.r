# nolint start: object_name_linter
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("C:\\Users\\Salma\\Documents\\Salma_Project_Data_All\\Kaib_paper\\C3_primary_model_julia_All_Mod\\FvCBestimation-master")
getwd()

# --- Step 1: Load saved results ---
scenario1_samples <- readRDS("scenario1_samples.rds")
scenario1_results <- readRDS("scenario1_results.rds")

# --- Step 2: Choose scenario, model, PPFD, and samples to plot ---
pct_range_str <- "10"  # Choose from "10", "20", "30", "50"
model_name <- "von Caemmerer (2000)"
ppfd_to_plot <- "500"
num_samples_to_plot <- 50

# --- Step 3: Extract simulation array ---
sim_results_key <- paste0("scenario1_", pct_range_str, "perc_results")
sim_list <- scenario1_results[[sim_results_key]]
if (is.null(sim_list)) stop("Scenario result not found!")
sim_array <- sim_list[[model_name]]
if (is.null(sim_array)) stop("Model not found in scenario results!")

# --- Step 4: Prepare data for plotting ---
ppfd_index <- which(dimnames(sim_array)$PPFD == ppfd_to_plot)
if (length(ppfd_index) == 0) stop("PPFD level not found!")
aci_mat <- sim_array[, ppfd_index, 1:num_samples_to_plot]
df <- as.data.frame(aci_mat)
df$Ci <- as.numeric(dimnames(sim_array)$Ci)
df_long <- pivot_longer(df, cols = -Ci, names_to = "Sample", values_to = "An")

# --- Step 5: Create plot ---
p <- ggplot(df_long, aes(x = Ci, y = An, group = Sample)) +
  geom_line(alpha = 0.3, color = "blue") +
  labs(
    title = paste("Scenario 1:", pct_range_str, "% parameter range"),
    subtitle = paste("Model:", model_name, "| PPFD:", ppfd_to_plot, "µmol m⁻² s⁻¹"),
    x = expression(paste("Intercellular CO"[2], " concentration (", mu, "mol mol"^-1, ")")),
    y = expression(paste("Net assimilation rate (", A[n], ", ", mu, "mol m"^-2, " s"^-1, ")"))
  ) +
  theme_minimal()

# --- Step 6: Save the plot ---
output_file <- paste0("ACi_curves_", model_name, "_PPFD", ppfd_to_plot, "_", pct_range_str, "perc.png")
ggsave(filename = output_file, plot = p, width = 8, height = 6, dpi = 300)

# Optionally print the plot to RStudio viewer as well
print(p)
