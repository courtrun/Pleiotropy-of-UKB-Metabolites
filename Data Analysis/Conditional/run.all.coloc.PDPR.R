set.seed(42)
library(ggplot2)
library(dplyr)

vars = c("none", "rs370014171_rs12924171", "bothalt", "three", "othertwo", "rs370014171", "rs8061221", "rs12924171")

traits = c("Ala", "Pyruvate", "Ile", "Val", "Leu")


data.table::rbindlist(lapply(vars, function(var) {
    data.table::rbindlist(lapply(traits, function(trait) {
        data.table::fread(paste0("conditional/PDPR/", var, ".", trait, ".glm.linear")) %>% mutate(trait = trait, condition=var)
    }), fill=T)
}), fill=T) %>% as.data.frame -> p
inputdata <- data.table::fread("../technical_variation/Residual_NMR_biomarkers.tsv", header=T)
library(matrixStats)
colSds(as.matrix(inputdata %>% select(-IID, -FID))) -> inputsds
sdcols <- data.frame(sdY = inputsds, trait = inputdata %>% select(-IID, -FID) %>% colnames)

inner_join(sdcols, p) -> sdp
sdp$type = "quant"
colnames(sdp) <- c("sdY", "trait", "chromosome", "position", "snp", "ref", "alt", "a1", "maf", "test", "n", "beta", "sebeta", "tstat", "p", "condition", "type")
library(coloc)
sdp %>% filter(condition == "none", !is.na(beta)) %>% mutate(varbeta = sebeta**2) %>% select(-condition) -> varp

# no clinical measures for PDPR
allp <- varp

data.table::rbindlist(lapply(unique(allp$trait), function(t1) {
    data.table::rbindlist(lapply(unique(allp$trait), function(t2) {
        d2 = as.list(allp %>% filter(trait == t2))
        d2$sdY = mean(d2$sdY)
        d2$type = unique(d2$type)
        d1 = as.list(allp %>% filter(trait == t1))
        d1$sdY = mean(d1$sdY)
        d1$type = unique(d1$type)
        as.data.frame(t(coloc.abf(dataset1 = d1, dataset2 = d2)$summary)) %>% mutate(t1=t1, t2=t2)
    }))
})) -> all.coloc

pdf("all.coloc.PDPR.pdf", useDingbats=F)
plt <- ggplot(all.coloc, aes(x=t1, y=t2, fill=PP.H4.abf)) + geom_tile() + theme_bw() + scale_fill_viridis_b()
print(plt);
dev.off()

write.table(all.coloc, "all_coloc_PDPR.tsv", quote=F, sep="\t", row.names=F, col.names=T)
