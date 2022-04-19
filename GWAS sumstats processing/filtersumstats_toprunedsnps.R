# Courtney Smith - BOLT-LMM - SNP pruning - Called in Snakefile_gwasprocessing script by rule pruned_snplist
# Goal of script: Filter across sumstats to just pruned snplist

library(dplyr)
### prune
bim <- data.table::fread("1000G_EUR_Phase3_plink/mymerge.bim", header=F)
colnames(bim) <- c("CHROM","SNP","CM","BP","REF","ALT")

args = commandArgs(trailingOnly=TRUE)

my.list <- list()
for (i in 1:(length(args)-2)){
filteredhits <- data.table::fread(commandArgs(TRUE)[i], fill=TRUE,header=T)
filteredhits <- filter(filteredhits,CHR!="CHR") %>% select(SNP,CHR,BP,ALLELE1,BETA,SE,P_BOLT_LMM)
filteredhits$CHR <- as.integer(filteredhits$CHR)
filteredhits$BP <- as.integer(filteredhits$BP)
filteredhits$P_BOLT_LMM <- as.numeric(filteredhits$P_BOLT_LMM)
my.list[[i]] <- filteredhits
}
filteredhits_combined <- bind_rows(my.list)
filteredhits_an <- inner_join(filteredhits_combined,bim,by=c("SNP","CHR"="CHROM","BP"))

fh <- as.data.frame(filteredhits_an %>% arrange(CHR,BP) %>% mutate(id=1:n(), prev=lag(CM, 1, 0), block=cumsum(CM - prev > 0.1)) %>% # gives ID to every variant and if it's more than 0.1 (was 0.01) centimorgans away from previous list then say it's a new region
     group_by(CHR, block) %>% filter(P_BOLT_LMM==min(P_BOLT_LMM))) # sample_n(1)
write.table(fh,commandArgs(TRUE)[length(args)-1], quote=F, sep="\t", row.names=F, col.names=T)

fh_s <- fh %>% select(SNP) %>% distinct()
write.table(fh_s,commandArgs(TRUE)[length(args)], quote=F, sep="\t", row.names=F, col.names=T)
