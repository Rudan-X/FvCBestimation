library(tidyverse)

# --- Table A (models) ---
tableA <- tribble(
  ~model_id, ~model_label,
  "1980FvCB",  "Farquhar et al. 1980",
  "1985Sharkey", "Sharkey et al. 1985",
  "1991Collatz", "Collatz et al. 1991",
  "1992aHarley", "Harley et al. 1992a",
  # "Harley1992b",  "Harley et al. 1992b",
  "2000vonCaemmerer",  "von Caemmerer et al. 2000",
  "2004Ethier", "Ethier & Livingston 2004",
  "2004Yin", "Yin et al. 2004",
  "2009Yin", "Yin et al. 2009",
  "2012Evans", "Evans & von Caemmerer 2012",
  "2012Tholen", "Tholen et al. 2012",
  "2016Moualeu",  "Moualeu et al. 2016",
  "2018Busch",  "Busch et al. 2018",
  "2020Busch",  "Busch et al. 2020"
)

# --- Table B (estimation studies) ---
tableB <- tribble(
  ~study_id, ~study_label, ~models_used, ~input_data,~method, ~software,
  "Harley1992b", "Harley et al. 1992b", "1980FvCB;1992aHarley","GE", "OLS", NA,
  # "Long2003", "Long & Bernacchi 2003", "FvCB1980","GE", "MLE", "SPSS",
  "Ethier2004", "Ethier & Livingston 2004", "2004Ethier","GE", "OLS", NA,
  "Dubois2007", "Dubois et al. 2007", "1980FvCB","GE", "OLS", "SAS",
  "Sharkey2007,16", "Sharkey et al. 2007, 2016", "2000vonCaemmerer","GE", "NLS", "Excel",
  "Su2009", "Su et al. 2009", "1980FvCB;2004Ethier","GE", "GA", "MATLAB",
  "Patrick2009", "Patrick et al. 2009", "2004Ethier","GE", "Hierarchical Bayesian", "WinBUGS",
  # "Yin2009", "Yin 2009", "Yin2004;Yin2009","GE + CF", "OLS variant", "SAS",
  "Gu2010", "Gu et al. 2010", "2000vonCaemmerer;2004Ethier","GE", "Hybrid search", "LeafWeb",
  "Yin2011", "Yin et al. 2011", "2009Yin","GE", "NLS", NA,
  "Sun2013", "Sun et al. 2013", "2004Ethier","GE", "NLS", NA,
  "Bellasio2015", "Bellasio et al. 2015", "2009Yin","GE + CF", "NLS", "Excel",
  # "Walker2015", "Walker & Ort 2015", "Ethier2004","GE", "Regression J", NA,
  "Plantecophys", "Plantecophys (Duursma et al. 2015)", "1991Collatz", "GE",  "NLS", "Rpackage",
  "Xiao2021", "Xiao et al. 2021", "2000vonCaemmerer","GE + CF", "Bayesion estimation", "Python",
  "msuRACiFit","msuRACiFit (Gregory et al. 2021)", "2000vonCaemmerer;2018Busch", "GE", "NLS", "Excel",
  "Photosynthesis","Photosynthesis (Stinziano1 et al. 2021)", "2000vonCaemmerer", "GE", "NLS", "Rpackage",
  "PhotoGEA","PhotoGEA (Lochocki et al. 2025)", "2000vonCaemmerer;2016Moualeu;2020Busch", "GE + CF", "NLS", "Rpackage",
  "PhotoTorch","PhotoTorch (Lei et al. 2025)", "2000vonCaemmerer;2004Ethier", "GE", "ADAM optimizer", "Python"

)
order_models <- c(
  "Farquhar et al. 1980",
  "Collatz et al. 1991",
  "Harley et al. 1992a",
  "von Caemmerer et al. 2000",
  "Ethier & Livingston 2004",
  "Yin et al. 2009",
  "Moualeu et al. 2016",
  "Busch et al. 2018",
  "Busch et al. 2020"
)


order_studies <- c(
  "Harley et al. 1992b",
  "Ethier & Livingston 2004",
  "Dubois et al. 2007",
  "Sharkey et al. 2007, 2016",
  "Su et al. 2009",
  "Patrick et al. 2009",
  "Gu et al. 2010",
  "Yin et al. 2011",
  "Sun et al. 2013",
  "Bellasio et al. 2015",
  "Plantecophys (Duursma et al. 2015)",
  "Xiao et al. 2021",
  "msuRACiFit (Gregory et al. 2021)",
  "Photosynthesis (Stinziano1 et al. 2021)",
  "PhotoGEA (Lochocki et al. 2025)",
  "PhotoTorch (Lei et al. 2025)"

)
# make sure the labels in your data EXACTLY match the strings above
crosswalk <- tableB %>%
  separate_rows(models_used, sep = ";") %>%
  left_join(tableA, by = c("models_used" = "model_id")) %>%
  select(study_id, study_label, method, software,
         model_id = models_used, model_label)
crosswalk <- crosswalk %>%
  mutate(
    model_label = factor(model_label, levels = rev(order_models)),
    study_label = factor(study_label, levels = order_studies)
  )
#
# ggplot(crosswalk, aes(x = model_label, y = study_label, fill = method)) +
#   geom_tile(color = "grey80") +
#   scale_x_discrete(limits = order_models, drop = FALSE) +  # enforce order
#   labs(x = "Model", y = "Estimation study", fill = "Method") +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# # --- Create crosswalk ---

# --- Show crosswalk
print(crosswalk)

# --- Optional: create a wide matrix (Study × Model) ---
matrix <- crosswalk %>%
  mutate(flag = 1) %>%
  pivot_wider(names_from = model_label, values_from = flag, values_fill = 0)

print(matrix)

# --- Visualization (heatmap-style) ---
ggplot(crosswalk, aes(x = study_label, y = model_label, fill = method)) +
  geom_tile(color = "grey80") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Estimation studies", y = "FvCB models", fill = "Method") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/Figures/model_vs_estimation.png",width = 12, height = 6)

