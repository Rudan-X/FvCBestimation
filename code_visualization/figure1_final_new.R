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
    FvCB80_Ac   [label='FvCB80:Ac' fillcolor='#FDE2E4']
    FvCB80_Aj   [label='FvCB80:Aj' fillcolor='#E6FFFB']
    FvCB80 [label='Farquhar (1980)']
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
    Busch18    [label='Busch (2018)']

    Harley92_Cc   [label='Harley92:Cc'  fillcolor='#E5D4FF']
    # Harley92_Aj   [label='Harley92:Aj'  fillcolor='#E6FFFB']
    Caemmerer00_Ap [label='Caemmerer00:Ap' fillcolor='#FFF1B6']
    Ethier04_Ac    [label='Ethier04:Ac'   fillcolor='#FDE2E4']
    Ethier04_Aj    [label='Ethier04:Aj'   fillcolor='#E6FFFB']
    Yin04_Aj       [label='Yin04:Aj'      fillcolor='#E6FFFB']
    Busch18_Aj     [label='Busch18:Aj'    fillcolor='#E6FFFB']
    Busch18_Ap     [label='Busch18:Ap'    fillcolor='#FFF1B6']
    Tholen12_gm    [label='Tholen12:Cc' fillcolor='#E5D4FF']


  }


  # FvCB components to other models
  FvCB80_Ac -> FvCB80 [color='black' style=solid] #label='inherits'
  FvCB80_Aj -> FvCB80 [color='black' style=solid] #label='inherits'

  FvCB80_Ac -> Harley92 [color='black' style=solid] #label='inherits'
  FvCB80_Aj -> Harley92 [color='black' style=solid]

  FvCB80_Ac -> Ethier04 [color='black' style=solid]
  FvCB80_Aj -> Ethier04 [color='black' style=solid]

  FvCB80_Ac -> Yin04 [color='black' style=solid]

  FvCB80_Ac -> Busch18 [color='black' style=solid]

  # FvCV Model to model
  FvCB80 -> Caemmerer00 [color='black' style=solid]
  FvCB80 -> Tholen12 [color='black' style=solid]

  # New model to new component
  Harley92 -> Harley92_Cc [label='gm' color='blue' style=dashed dir=back]
  # Harley92 -> Harley92_Aj [label='gm' color='blue' style=dashed]

  Caemmerer00 -> Caemmerer00_Ap [label='TPU from Sharkey' color='blue' style=dashed  dir=back]

  Caemmerer00 -> Dubois07 [label='min(Ac,Aj,Ap)']

  Ethier04 -> Ethier04_Ac [label='Quadratic form' color='blue' style=dashed  dir=back]
  Ethier04 -> Ethier04_Aj [label='Quadratic form' color='blue' style=dashed  dir=back]
  Ethier04 -> Xiao20 [color='black' style=solid]

  Yin04 -> Yin04_Aj [label='Alternative ETC' color='blue' style=dashed  dir=back]

  # Dubois07 -> Dubois07_gm [label='gm & leakiness' color='blue' style=dashed]

  Tholen12 -> Tholen12_gm [label='New gm definition (gch,gwp)' color='blue' style=dashed  dir=back]

  Busch18 -> Busch18_Aj [label='N assimilation' color='blue' style=dashed  dir=back]
  Busch18 -> Busch18_Ap [label='N assimilation' color='blue' style=dashed  dir=back]


}
")

svg_txt <- export_svg(g)

# 3. Convert to PNG (high resolution)
rsvg_png(charToRaw(svg_txt),
         file = "results/Figures/Fig1.Model_diagram.png",
         width = 2200, height = 1500)

#############################################################
load(file="results/Rdata/fig1_simulated_curves_nonoise.RData")


curves_name <- c("FvCB80", "Harley92","vonC00",
                 "Ethier04", "Yin04", "Dubois07",
                 "Tholen12", "Busch18", "Xiao21" )
names(simulated_curves) <- curves_name

df <- map_dfr(simulated_curves, ~as.data.frame(.x), .id = "Model")
df$PPFD <- paste0("PPFD: ",df$PPFD)
df$PPFD <- factor(df$PPFD, levels=c("PPFD: 50","PPFD: 150","PPFD: 300" ,"PPFD: 500", "PPFD: 1200" ,"PPFD: 1800"))
df$Model <- factor(df$Model,levels=curves_name)

ggplot(df, aes(Ci, An, color = Model, shape=Model)) + #
  # geom_line(linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 1.5 ) +
  facet_wrap(~PPFD, ncol = 3, dir = "h") +
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
