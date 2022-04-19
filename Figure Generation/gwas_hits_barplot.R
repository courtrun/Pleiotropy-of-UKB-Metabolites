# Courtney Smith - Nightingale Analysis - BOLT-LMM - GWAS hits barplot
# Goal of script: Plot the number of GWAS hits for each met

library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)

# Set threshold of sig and metabolites to include
sig_thres = 5e-8
keep_rigid <- c("Gly","Ala","Leu","Val","Ile","Pyruvate","Tyr","His","Phe","Gln","Acetone","bOHbutyrate","Acetoacetate","Glucose","Lactate","Citrate")

# coloring
cat <- data.table::fread("biomarkers_pathwayorder.csv") %>% select(Name_in_TSV,Color_Group,Biomarker_name)
cat <- as.data.frame(cat)
# Add BCAA category
BCAA_list <- c("Leu","Ile","Val")
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[1]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[2]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[3]] <- "BCAA"
m <- filter(cat,Name_in_TSV %in% keep_rigid) %>% mutate(met_color_label=ifelse(Color_Group=="Amino_Acid","#ffaa00",ifelse(Color_Group=="BCAA","#00aaff",ifelse(Color_Group=="Glycolysis","#ff0000","#aa00ff"))))

# Number of clumped hits in each met
mets <- list()
for (met in keep_rigid) {
num_hits <- length(unique(data.table::fread(paste0("metabolomeFirstPass.",met,".hits"))$SNP))
mets[[met]] <- data.frame(Metabolite=met,num_hits=num_hits)
}
clump_met_hits <- bind_rows(mets) %>% left_join(m,by=c("Metabolite"="Name_in_TSV")) %>% arrange(Color_Group,desc(num_hits)) %>% mutate(order=head(letters,nrow(.)))

ggplot(clump_met_hits,aes(x=fct_reorder(Biomarker_name,order),y=num_hits,fill=Color_Group))+geom_bar(stat="identity",position="dodge")+
  labs(y="Number of GWAS hits",x="Metabolite",fill="Metabolite Group")+
    scale_fill_manual(values=unique(clump_met_hits$met_color_label))+
    theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text=element_text(size=16))
ggsave("clump_GWAS_hits.pdf")
