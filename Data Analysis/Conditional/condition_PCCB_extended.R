library(ggplot2)
library(dplyr)

#vars = c("none", "both", "bothalt", "three", "othertwo", "rs12325419", "rs151047766", "rs34329336", "rs370014171", "rs4985532", "rs554468395", "rs56375022", "rs6499295", "rs6499327", "rs8061221", "rs12924171")
vars = c("none", "rs645040", "rs61791721", "rs1247764225", "both", "rs687339")

traits = c("Ala", "Pyruvate", "Gly", "Ile", "Val", "Leu", "Total_FA", "Total_FC", "Total_TG", "LDL_C", "HDL_C", "CAD")

data.table::rbindlist(lapply(vars, function(var) {
    data.table::rbindlist(lapply(traits, function(trait) {
        if (trait == "CAD") {
        data.table::fread(paste0("conditional/PCCBextended/", var, ".", trait, ".glm.logistic.hybrid")) %>% mutate(trait = trait, condition=var)
        } else {
        data.table::fread(paste0("conditional/PCCBextended/", var, ".", trait, ".glm.linear")) %>% mutate(trait = trait, condition=var)
        }
    }), fill=T)
}), fill=T) %>% as.data.frame -> p
        
pdf("condition_PCCBextended.pdf", useDingbats=F, height=2*(length(traits) + 1), width=2*(length(vars)+1))
g <- ggplot(p %>% filter(A1_FREQ > 1e-2, A1_FREQ < 0.99), aes(x=POS, y=-log10(P))) + geom_point() + theme_bw() + facet_grid(trait ~ condition, scale="free")
print(g)
dev.off()
