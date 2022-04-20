# Courtney Smith - Nightingale Analysis - BOLT-LMM - Metabolite GWAS hit gene distance distribution
# Goal of script: Distance of each variant to assigned gene by variant type

library(dplyr)
library(ggplot2)
library(ggridges)

s <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")

ggplot(s, aes(x = Dist_Assigned_Gene_Min, y = Variant_Classification, fill = Variant_Classification)) +
  geom_density_ridges() + labs(x="Distance to Gene",y="Variant Classification")+
  theme_ridges() +  theme(legend.position = "none")
ggsave("rigidline_distminassignedgene_byvt.png")
