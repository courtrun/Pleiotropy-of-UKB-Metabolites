
library(ggplot2)
library(dplyr)

#vars = c("none", "both", "bothalt", "three", "othertwo", "rs12325419", "rs151047766", "rs34329336", "rs370014171", "rs4985532", "rs554468395", "rs56375022", "rs6499295", "rs6499327", "rs8061221", "rs12924171")
vars = c("none", "rs370014171_rs12924171", "bothalt", "three", "othertwo", "rs370014171", "rs8061221", "rs12924171")

traits = c("Ala", "Pyruvate", "Ile", "Val", "Leu")
#traits = c("Ala", "Pyruvate", "Ile", "Val", "Leu", "Lactate", "Gln", "bOHbutyrate")

data.table::rbindlist(lapply(vars, function(var) {
    data.table::rbindlist(lapply(traits, function(trait) {
        data.table::fread(paste0("conditional/PDPR/", var, ".", trait, ".glm.linear")) %>% mutate(trait = trait, condition=var)
    }))
})) %>% as.data.frame -> p
        
pdf("condition_PDPR.pdf", useDingbats=F, height=2*(length(traits) + 1), width=2*(length(vars)+1))
g <- ggplot(p %>% filter(A1_FREQ > 1e-2, A1_FREQ < 0.99), aes(x=POS, y=-log10(P))) + geom_point() + theme_bw() + facet_grid(trait ~ condition, scale="free")
print(g)
dev.off()
