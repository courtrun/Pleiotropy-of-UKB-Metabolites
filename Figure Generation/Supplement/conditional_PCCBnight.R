# Courtney Smith - Nightingale Analysis - BOLT-LMM - Conditional analysis - PCCB
# Goal of script: Plot conditional analysis results for PCCB with nightingale metabolites

library(ggplot2)
library(dplyr)

vars = c("none", "rs645040", "rs61791721")
traits = c("Ala", "Pyruvate", "Gly", "Ile", "Val", "Leu", "Total_FA", "Total_TG") # "Ala", "Pyruvate", "Gly", "Ile", "Val", "Leu", "Total_FA", "Total_FC", "Total_TG", "LDL_C", "HDL_C", "CAD"

data.table::rbindlist(lapply(vars, function(var) {
    data.table::rbindlist(lapply(traits, function(trait) {
        if (trait == "CAD") {
        data.table::fread(paste0("conditional/PCCBextended/", var, ".", trait, ".glm.logistic.hybrid")) %>% mutate(trait = trait, condition=var)
        } else {
        data.table::fread(paste0("conditional/PCCBextended/", var, ".", trait, ".glm.linear")) %>% mutate(trait = trait, condition=var)
        }
    }), fill=T)
}), fill=T) %>% as.data.frame -> p

# Update trait names
keep_rigid <- data.frame(Met=c("bOHbutyrate","Acetoacetate","Acetone","Citrate","Lactate","Pyruvate","Phe","Tyr","Glucose","Ala","Leu","Ile","Val","Gly","Gln","His", "Total_FA", "Total_TG"))
m <- data.frame(Biomarker_name=c("3-Hydroxybutyrate","Acetoacetate","Acetone","Citrate","Lactate","Pyruvate","Phenylalanine","Tyrosine","Glucose","Alanine","Leucine","Isoleucine","Valine","Glycine","Glutamine","Histidine", "Total Fatty Acids", "Total Triglycerides"))
m <- bind_cols(m,keep_rigid)

# Update condition/variant names
v <- data.frame(conditions=c("Marginal","rs645040", "rs61791721"))
v <- bind_cols(v,distinct(p %>% select(condition)))

p <- left_join(p,m,by=c("trait"="Met"))
p <- left_join(p,v,by=c("condition"))
p <- filter(p,conditions!= "rs645040")
p$conditions = factor(p$conditions, levels=c("Marginal","rs61791721"))
p$Biomarker_name = factor(p$Biomarker_name, levels=c("Alanine","Leucine","Isoleucine","Valine","Pyruvate","Glycine","Total Fatty Acids","Total Triglycerides"))

g <- ggplot(p %>% filter(A1_FREQ > 1e-2, A1_FREQ < 0.99), aes(x=POS, y=-log10(P))) +
      geom_point() + theme_bw() + labs(x=paste0("Position on Chromosome ",unique(p$`#CHROM`)))+
      theme(axis.text.x = element_text(angle = 70, hjust=1))+
      facet_grid(Biomarker_name ~ conditions, scale="free")

pdf("condition_PCCBnight.pdf", useDingbats=F, height=2*(length(traits) + 1), width=2*(length(vars)+1))
print(g)
dev.off()
