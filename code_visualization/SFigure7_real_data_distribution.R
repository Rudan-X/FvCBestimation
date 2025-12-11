library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(ggpubr)

setwd("C:/Users/Rudan/Documents/GitHub/FvCBestimation/")

ACi <- read.csv("data/2022_ACi_rawData_barley.csv")
AQ <- read.csv("data/2022_AQcurves_rawData_barley.csv")
ref <- read.csv("data/2022_barley_reference.csv")
AQ <- AQ[AQ$Flag_removal!="x",]

ACi$ID <- paste0(ACi$Plot,"_",ACi$Repeat)
AQ$ID <- paste0(AQ$Plot,"_",ACi$Repeat)
ACi$genotype <- ref$Accession[match(ACi$Plot,ref$PlotID)]
AQ$genotype <- ref$Accession[match(AQ$Plot,ref$PlotID)]

rmind <- which(ACi$genotype=="B1K-14-17" & ACi$Photo < (-5))
ACi <- ACi[-rmind,]

genotype_id <- intersect(ACi$genotype, AQ$genotype)
genotype_id <- genotype_id[1:15]
# genotype_id <- c("B1K-03-17","B1K-12-10","B1K-25-08","B1K-23-19","B1K-26-13",
#                  "B1K-27-19","B1K-28-13","B1K-30-13","B1K-36-15","B1K-45-05")


ACi <- ACi[ACi$genotype%in%genotype_id,c("Ci","Photo","PPFD","genotype","ID", "Plot")]
AQ <- AQ[AQ$genotype%in%genotype_id,c("Ci","A","PAR","genotype","ID", "Plot")]


ACi <- ACi %>%
  group_by(genotype) %>%
  mutate(
    ID_new = match(ID, unique(ID)),
    ID_new = as.character(ID_new)
  ) %>%
  ungroup()

AQ <- AQ %>%
  group_by(genotype) %>%
  mutate(
    ID_new = match(ID, unique(ID)),
    ID_new = as.character(ID_new)
  ) %>%
  ungroup()

g1 <- ggplot(ACi, aes(x = Ci, y = Photo,  color = ID_new)) +
  geom_point() +
  facet_wrap(~ genotype, ncol=5) +
  theme_bw() +
  labs(x = expression(C[i]~(µmol~mol^{-1})),
       y = expression(A[n]~(µmol~m^{-2}~s^{-1})),
       color = "Plant") +
  theme(legend.position = "none")

g1

g2 <- ggplot(AQ, aes(x = PAR, y = A, color = ID_new)) +
  geom_point() +
  facet_wrap(~ genotype, ncol=5) +
  theme_bw() +
  labs(x = expression(PPFD~(µmol~m^{-2}~s^{-1})),
       y = expression(A[n]~(µmol~m^{-2}~s^{-1})),
       color = "Plant") +
  theme(legend.position = "bottom")

g2

ggarrange(g1,g2,labels = c('a', 'b'), ncol=1, heights = c(0.85,1))

ggsave(filename = paste0("results/Figures/SFig7_real_data_distribution.png"),width = 7, height = 9)



