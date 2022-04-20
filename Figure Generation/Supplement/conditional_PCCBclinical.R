Courtney Smith - Nightingale Analysis - BOLT-LMM - Conditional analysis - PCCB
# Goal of script: Plot conditional analysis results for PCCB with ukbb clinical measures

library(ggplot2)
library(dplyr)

vars = c("none", "rs61791721", "rs1247764225", "both") # rs645040

traits = c("HDL_cholesterol", "Triglycerides","LDL_direct", "softNoSR", "hardNoSR", "soft", "hard") # "LDL_direct", "Cholesterol", "Urea", "HDL_cholesterol", "SHBG", "Testosterone", "CAD", "Triglycerides",, "softNoSR", "hardNoSR", "soft", "hard"

data.table::rbindlist(lapply(vars, function(var) {
    data.table::rbindlist(lapply(traits, function(trait) {
        if (trait %in% c("softNoSR", "hardNoSR", "soft", "hard")) {
        data.table::fread(paste0("conditional/PCCBclinical/", var, ".", trait, ".glm.logistic.hybrid")) %>% mutate(trait = trait, condition=var)
        } else {
        data.table::fread(paste0("conditional/PCCBclinical/", var, ".", trait, ".glm.linear")) %>% mutate(trait = trait, condition=var)
        }
    }), fill=T)
}), fill=T) %>% as.data.frame -> p

# Update trait names
mets <- data.frame(Met=c("HDL_cholesterol", "Triglycerides","LDL_direct","softNoSR", "hardNoSR", "soft", "hard"))
m <- data.frame(Biomarker_name=c("HDL_C","Triglycerides","LDL_C","softNoSR", "hardNoSR", "soft", "CAD"))
m <- bind_cols(m,mets)

# Update condition/variant names
v <- data.frame(conditions=c("Marginal","rs61791721", "rs1247764225", "rs61791721 and rs1247764225"))
v <- bind_cols(v,distinct(p %>% select(condition)))

p <- left_join(p,m,by=c("trait"="Met"))
p <- left_join(p,v,by=c("condition"))
p <- filter(p,Biomarker_name %in% c("HDL_C", "Triglycerides","LDL_C","CAD"))
p$conditions = factor(p$conditions, levels=c("Marginal","rs61791721", "rs1247764225", "rs61791721 and rs1247764225"))
p$Biomarker_name = factor(p$Biomarker_name, levels=c("HDL_C","LDL_C", "Triglycerides","CAD"))

g <- ggplot(p %>% filter(A1_FREQ > 1e-2, A1_FREQ < 0.99), aes(x=POS, y=-log10(P))) +
      geom_point() + theme_bw() + labs(x=paste0("Position on Chromosome ",unique(p$`#CHROM`)))+
      theme(axis.text.x = element_text(angle = 70, hjust=1))+
      facet_grid(Biomarker_name ~ conditions, scale="free")

pdf("condition_PCCBclinical_CAD.pdf", useDingbats=F, height=2*(length(traits) + 1), width=2*(length(vars)+1))
print(g)
dev.off()
