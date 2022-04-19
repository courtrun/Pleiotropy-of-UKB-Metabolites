# Modified by Courtney Smith from Nasa's script - Called in Snakefile_gwasprocessing by rule filter_plink
# Goal of script: Filters by MAF and INFO

library(dplyr)
sumstats <- data.table::fread(commandArgs(TRUE)[1], fill=TRUE,header=T)

chromosome.order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "XY", "Y", "M")

temp <- sumstats %>% filter(A1FREQ >= 0.01, A1FREQ <= 0.99, INFO > 0.7, !is.na(P_BOLT_LMM))
mfs <- temp %>% arrange(factor(as.character(CHR), levels=chromosome.order), BP)

write.table(mfs, gzfile(commandArgs(TRUE)[2]), quote=F, sep="\t", row.names=F, col.names=T)
