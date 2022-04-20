# Courtney Smith - Nightingale Analysis - BOLT-LMM - Conditional analysis - SLC36A2
# Goal of script: Plot conditional analysis results for SLC36A2

library(ggplot2)
library(dplyr)

vars = c("none", "rs77010315")

traits = c("Ala", "Pyruvate", "Ile", "Val", "Leu")

data.table::rbindlist(lapply(vars, function(var) {
    data.table::rbindlist(lapply(traits, function(trait) {
        data.table::fread(paste0("conditional/SLC36A2/", var, ".", trait, ".glm.linear")) %>% mutate(trait = trait, condition=var)
    }))
})) %>% as.data.frame -> p

# Update trait names
keep_rigid <- data.frame(Met=c("bOHbutyrate","Acetoacetate","Acetone","Citrate","Lactate","Pyruvate","Phe","Tyr","Glucose","Ala","Leu","Ile","Val","Gly","Gln","His"))
m <- data.frame(Biomarker_name=c("3-Hydroxybutyrate","Acetoacetate","Acetone","Citrate","Lactate","Pyruvate","Phenylalanine","Tyrosine","Glucose","Alanine","Leucine","Isoleucine","Valine","Glycine","Glutamine","Histidine"))
m <- bind_cols(m,keep_rigid)

# Update condition/variant names
v <- data.frame(conditions=c("Marginal","rs77010315"))
v <- bind_cols(v,distinct(p %>% select(condition)))

p <- left_join(p,m,by=c("trait"="Met"))
p <- left_join(p,v,by=c("condition"))
p$conditions = factor(p$conditions, levels=c("Marginal","rs77010315"))

g <- ggplot(p %>% filter(A1_FREQ > 1e-2, A1_FREQ < 0.99), aes(x=POS, y=-log10(P))) +
      geom_point() + theme_bw() + labs(x=paste0("Position on Chromosome ",unique(p$`#CHROM`)))+
      theme(axis.text.x = element_text(angle = 70, hjust=1))+
      facet_grid(Biomarker_name ~ conditions, scale="free")

pdf("condition_SLC36A2.pdf", useDingbats=F, height=2*(length(traits) + 1), width=2*(length(vars)+1))
print(g)
dev.off()
