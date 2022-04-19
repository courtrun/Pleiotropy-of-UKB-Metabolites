# Courtney Smith - BOLT-LMM - SNP by Metabolite matrix - Called in Snakefile_gwasprocessing by rule make_matrix
# Goal of script: Make matrix of pruned hits by rigid metabolites

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

snps <- data.table::fread(commandArgs(TRUE)[1], fill=TRUE,header=T)
snps <- unique(snps$SNP)

# Combine hits for all rigid mets
my.list <- list()
for (i in 2:(length(args)-2)){
met_name=gsub("_.*","",gsub(".*filtered/","",commandArgs(TRUE)[i]))
my.list[[i]] <- data.table::fread(commandArgs(TRUE)[i],fill=TRUE)%>% mutate(TraitName=met_name) %>%
 mutate(P_BOLT_LMM=as.numeric(P_BOLT_LMM))%>% select(SNP,CHR,BP,ALLELE1,BETA,SE,P_BOLT_LMM,TraitName)
}
rigidmethits <- bind_rows(my.list)
met_long <- filter(rigidmethits,SNP %in% snps)

# Filter to only keep duplicate SNPs' row with lowest p-value
library(data.table)
met_long <- as.data.table(met_long %>% group_by(SNP,TraitName) %>% filter(P_BOLT_LMM == min(P_BOLT_LMM))%>% filter(1:n() == 1)) %>% rename(CHROM=CHR)
met_long$P_BOLT_LMM <- as.numeric(met_long$P_BOLT_LMM)
write.table(met_long,commandArgs(TRUE)[length(args)-1], quote=F, sep="\t", row.names=F, col.names=T)

# Long to wide form
met_wide <- dcast(as.data.table(met_long), SNP + CHROM + BP ~ paste0("_", TraitName), fun=mean,value.var = c("BETA","SE","P_BOLT_LMM"),sep='')
write.table(met_wide,commandArgs(TRUE)[length(args)], quote=F, sep="\t", row.names=F, col.names=T)
