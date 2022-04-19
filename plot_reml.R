# Courtney Smith - Nightingale Analysis - BOLT-LMM - Local Gen Cor barplot
# Goal of script: Plot bars for local gen cor

library(dplyr)
library(ggplot2)

gc <- data.table::fread("bolt_reml_Ala_Gln_estimates.tsv")
gc <- as.data.frame(gc)
gc["V1"][gc["V1"]=="CUSTOM_GLYCOLYSIS_GLUCONEOGENESIS_TCA"] <- "Glycolysis, Gluconeogenesis\nand TCA Genes"
gc["V1"][gc["V1"]=="CUSTOM_AA"] <- "Other Amino Acid Genes"
gc["V1"][gc["V1"]=="CUSTOM_BCAA"] <- "BCAA Genes"
gc["V1"][gc["V1"]=="CUSTOM_UREAMET"] <- "Urea Cycle Genes"
gc["V1"][gc["V1"]=="LIPID"] <- "Lipid Genes"
gc["V1"][gc["V1"]=="CUSTOM_KETONE"] <- "Ketone Body Genes"
gc["V1"][gc["V1"]=="CUSTOM_TF"] <- "Metabolite Associated\nTF Genes"
gc["V1"][gc["V1"]=="CUSTOM_TRANSPORTERS"] <- "Metabolite Associated\nTransporter Genes"
gc["V1"][gc["V1"]=="Coding"] <- "All Regions"
gc$V1 <- as.factor(gc$V1)
gc <- filter(gc,V1!="Lipid Genes")
gc <- filter(gc,V1!="CUSTOM_RIGID")
gc <- filter(gc,V1!="Ketone Body Genes")
gc <- gc %>% mutate(ordering=ifelse(grepl("All Regions",V1),100,ifelse(grepl("Metabolite",V1),110,V2)))
color_list <- c("All Regions"="gray","Other Amino Acid Genes"="#ffaa00","BCAA Genes"="#00aaff","Glycolysis, Gluconeogenesis\nand TCA Genes"="#ff0000",
      "Metabolite Associated\nTF Genes"="yellow4","Metabolite Associated\nTransporter Genes"="#ea70ff","Urea Cycle Genes"="#00b881") # "#a5751d","pink","Ketone Body Genes"="#aa00ff",

ggplot(gc,aes(x=reorder(V1,ordering), y=V2, ymin=V2-V3, ymax=V2+V3, fill=V1)) +
  geom_bar(stat="identity") + labs(y="Local Genetic Correlation for Alanine and Glutamine",x="Genome Region")+
  geom_errorbar() + theme_classic() + #scale_fill_brewer(palette="Dark2") +
   coord_flip() + scale_fill_manual(values=color_list)+
   geom_vline(xintercept=4.5,linetype="longdash")+
  guides(fill = "none")+ theme(axis.text.y = element_text(size=20),axis.text.x = element_text(size=24),
    axis.title.x = element_text(size=24),axis.title.y = element_text(size=24))
ggsave("bolt_reml_Ala_Gln_localgencor.pdf",width=20,height=10)
