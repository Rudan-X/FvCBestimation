library(tidyverse)
library(ggrepel)

tableA <- tribble(
  ~model_id, ~model_label, ~contribution, ~contri_detail,
  "1980FvCB","Farquhar et al. 1980","An definition","min(Ac, Aj)",
  "1985Sharkey","Sharkey et al. 1985","TPU limitation","min(Ac, Aj, Ap)",
  "1991Collatz","Collatz et al. 1991","An definition","Quadratic smooth of Ac, Aj",
  "1992aHarley","Harley et al. 1992a","Mesophyll conductance","Cc=Ci-A/gm",
  "2000vonCaemmerer","von Caemmerer et al. 2000","An definition","Review book",
  "2004Ethier","Ethier & Livingston 2004","An definition","Quadratic formulation of A(Ci,gm)",
  "2004Yin","Yin et al. 2004","RuBP-regeneration","Non-linear ETC",
  "2009Yin","Yin et al. 2009","RuBP-regeneration;Mesophyll conductance", "J(PPFD,Phi2),gm estimation",
  "2012Evans","Evans & von Caemmerer 2012","Mesophyll conductance","gm(g_liq, g_mem)",
  "2012Tholen","Tholen et al. 2012","Mesophyll conductance","rm(r_wp,r_ch)",
  "2016Moualeu","Moualeu et al. 2016","Mesophyll conductance","gm calculated using CF",
  "2018Busch","Busch et al. 2018","RuBP-regeneration;TPU limitation","alpha_S",
  "2020Busch","Busch et al. 2020","RuBP-regeneration;TPU limitation","alpha_S,alpha_T"
)


df <- tableA %>%
  mutate(year = as.numeric(str_extract(model_id, "^\\d{4}"))) %>%
  arrange(year) %>%
  mutate(model_label = factor(model_label, levels = model_label))

ggplot(df, aes(x = 0, y = year, color = contribution)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = model_label),
                  nudge_x = 0.25, size = 4, segment.color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.6, color = "grey80") +
  scale_y_reverse(breaks = unique(df$year)) +   # <-- reverse order
  labs(x = NULL, y = NULL, color = "Contribution") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave("results/Figures/FvCBmodels.png",width = 12, height = 6)

###############################
df2 <- df %>%
  mutate(contribution = if_else(contribution == "", "(unspecified)", contribution))

# order x by year
x_order <- df2 %>% arrange(year) %>% pull(model_label) %>% unique()

ggplot(df2, aes(x = contribution,
                y = factor(model_label, levels = x_order), fill = contribution)) +
  geom_tile(width = 0.95, height = 0.9, color = "white") +
  geom_text(aes(label = ifelse(contri_detail == "", "\u2013", contri_detail)),
            size = 3) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  guides(fill = "none") +
  labs(x = NULL, y = NULL,
       title = "FvCB-related contributions by study",
       subtitle = "Text inside tiles shows the specific formulation/detail") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA),
        legend.position = "bottom")


