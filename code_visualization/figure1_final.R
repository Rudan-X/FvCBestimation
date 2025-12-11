library(DiagrammeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(ggpubr)


library(DiagrammeRsvg)  # converts grViz -> SVG
library(rsvg)            # converts SVG -> PNG/PDF


g <- grViz("
digraph {
  graph [layout = dot, rankdir = LR]

  node [shape=box style=filled fontname=Helvetica]

  subgraph cluster_level1 { rank=same; style=invis;
    FvCB98_Ac   [label='FvCB98:Ac' fillcolor='#FDE2E4']
    FvCB98_Aj   [label='FvCB98:Aj' fillcolor='#E6FFFB']
    FvCB98 [label='Farquhar (1980)']
  }

  ## Level 2: Models (gray)
  subgraph cluster_level2 {
    Harley92  [label='Harley (1992)']
    Caemmerer00[label='von Caemmerer (2000)']
    Ethier04   [label='Ethier &amp; Livingston (2004)']
    Yin04      [label='Yin (2004)']
    Dubois07    [label='Dubois (2007)']
    Tholen12   [label='Tholen (2012)']
    Xiao20    [label='Xiao (2021)']
    Busch18    [label='Busch (2017)']

    Harley92_Cc   [label='Harley92:Cc'  fillcolor='white']
    # Harley92_Aj   [label='Harley92:Aj'  fillcolor='#E6FFFB']
    Caemmerer00_Ap [label='Caemmerer00:Ap' fillcolor='#FFF1B6']
    Ethier04_Ac    [label='Ethier04:Ac'   fillcolor='#FDE2E4']
    Ethier04_Aj    [label='Ethier04:Aj'   fillcolor='#E6FFFB']
    Yin04_Aj       [label='Yin04:Aj'      fillcolor='#E6FFFB']
    Busch18_Aj     [label='Busch17:Aj'    fillcolor='#E6FFFB']
    Busch18_Ap     [label='Busch17:Ap'    fillcolor='#FFF1B6']
    Tholen12_gm    [label='Tholen12:gch,gwp' fillcolor='#E5D4FF']


  }


  # FvCB components to other models
  FvCB98_Ac -> FvCB98 [color='black' style=solid] #label='inherits'
  FvCB98_Aj -> FvCB98 [color='black' style=solid] #label='inherits'

  FvCB98_Ac -> Harley92 [color='black' style=solid] #label='inherits'
  FvCB98_Aj -> Harley92 [color='black' style=solid]

  FvCB98_Ac -> Ethier04 [color='black' style=solid]
  FvCB98_Aj -> Ethier04 [color='black' style=solid]

  FvCB98_Ac -> Yin04 [color='black' style=solid]

  FvCB98_Ac -> Busch18 [color='black' style=solid]

  # FvCV Model to model
  FvCB98 -> Caemmerer00 [color='black' style=solid]
  FvCB98 -> Tholen12 [color='black' style=solid]

  # New model to new component
  Harley92 -> Harley92_Cc [label='gm' color='blue' style=dashed]
  # Harley92 -> Harley92_Aj [label='gm' color='blue' style=dashed]

  Caemmerer00 -> Caemmerer00_Ap [label='TPU from Sharkey' color='blue' style=dashed]

  Caemmerer00 -> Dubois07 [label='min(Ac,Aj,Ap)' color='blue' style=dashed]

  Ethier04 -> Ethier04_Ac [label='Quadratic form' color='blue' style=dashed]
  Ethier04 -> Ethier04_Aj [label='Quadratic form' color='blue' style=dashed]
  Ethier04 -> Xiao20 [color='black' style=solid]

  Yin04 -> Yin04_Aj [label='Alternative ETC' color='blue' style=dashed]

  # Dubois07 -> Dubois07_gm [label='gm & leakiness' color='blue' style=dashed]

  Tholen12 -> Tholen12_gm [label='New gm definition' color='blue' style=dashed]

  Busch18 -> Busch18_Aj [label='N assimilation' color='blue' style=dashed]
  Busch18 -> Busch18_Ap [color='blue' style=dashed]


}
")

svg_txt <- export_svg(g)

# 3. Convert to PNG (high resolution)
rsvg_png(charToRaw(svg_txt),
         file = "results/Figures/Fig1.Model_diagram.png",
         width = 2200, height = 1500)

#############################################################
load(file="results/Rdata/fig1_simulated_curves_nonoise.RData")


curves_name <- c("FvCB80", "Harley92","vonCaemmerer00",
                 "Ethier04", "Yin04", "Dubois07",
                 "Tholen12", "Busch17", "Xiao21" )
names(simulated_curves) <- curves_name

df <- map_dfr(simulated_curves, ~as.data.frame(.x), .id = "Model")
df$PPFD <- paste0("PPFD: ",df$PPFD)
df$PPFD <- factor(df$PPFD, levels=c("PPFD: 50","PPFD: 150","PPFD: 300" ,"PPFD: 500", "PPFD: 1200" ,"PPFD: 1800"))
df$Model <- factor(df$Model,levels=curves_name)

ggplot(df, aes(Ci, An, color = Model, shape=Model)) + #
  # geom_line(linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 1.5 ) +
  facet_wrap(~PPFD, ncol = 3, dir = "v") +
  labs(x = expression(C[i]~"(µmol mol"^{-1}*")"),
       y = expression(An~"(µmol m"^{-2}~s^{-1}*")"),
       color = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")+
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 18, 7, 10))  # c

ggsave("results/Figures/Fig1.simulation.png",width = 7, height = 6)

library(cowplot)
library(magick)

# Read as raster images
p1 <- ggdraw() + draw_image("results/Figures/Fig1.Model_diagram.png")
p2 <- ggdraw() + draw_image("results/Figures/Fig1.simulation.png")

# Combine horizontally
combined <- plot_grid(p1, p2, ncol = 1, align = "v", labels = c("a","b"), rel_heights = c(0.8,1))
combined


# Save combined figure
ggsave("results/Figures/Fig1_combined.png", combined, width = 7, height = 9.5, dpi = 300)
