# Courtney Smith - Nightingale Analysis - BOLT-LMM - Discordant Variant Analysis
# Goal of script: Data table for discordant and concordant pathway-relevant enzyme variant annotation

library(dplyr)

d <- data.table::fread("disvar_annotation.tsv")
dfilt <- d %>% filter(Assigned_Gene_Type=="Pathway_Relevant_Enzyme") %>% arrange(desc(Assigned_Gene)) %>% select(SNP,Assigned_Gene,Assigned_Gene_Type,Metabolite.1,Metabolite.2)
dfilt <- as.data.frame(dfilt %>% mutate(Relative_Location="Between"))
dfilt["Relative_Location"][dfilt["Assigned_Gene"]=="PFKP"] <- "Not_Between"

# Annotate the discordant variants

# Annotate concordant variants
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")
an <- filter(an,Variant_Classification=="Concordant")

library("ashr")
library("dplyr")

gencor <- data.table::fread("combined_gen_cor_sumstats.txt",fill=TRUE)
gencor <- na.omit(gencor) # get rid of pairs where one or more files did not exist
# Reformat trait columns to just be metabolite names
gencor$p1 <- gsub(".*Pass.","",gencor$p1)
gencor$p1 <- gsub(".txt.*","",gencor$p1)
gencor$p2 <- gsub(".*Pass.","",gencor$p2)
gencor$p2 <- gsub(".txt.*","",gencor$p2)
gencor <- gencor %>% filter(p1!=p2) # Drop pairs of metabolites with themselves
gencor <- gencor[!duplicated(t(apply(gencor %>% select(p1,p2), 1, sort))),] # get rid of repeat pairs

# filter to pairs where both mets are rigid mets
full <- data.table::fread("mGWAS_rigidmetabolites_sig_pruned.tsv")
rigidmets_list <- gsub("BETA_","",colnames(full %>% select(matches("BETA_"))))
gc <- filter(gencor,p1 %in% rigidmets_list & p2 %in% rigidmets_list)
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")

# filter to pairs sig dif than 0
filt <- ashr::ash(gc$rg,gc$se)
lfsr <- as.data.frame(filt$result) %>% select(lfsr)
gencor_lfsr <- cbind(gc,lfsr)
gencor_lfsr <- filter(gencor_lfsr,lfsr<0.005) # 0.005 because that's less than 5 in 1000 and theres's about 120 correlations

# Identify variants sig in at least one rigid metabolite 5e-8 and 1e-4 in another
pvals <- full %>% select(matches("P_BOLT_LMM_"))
sig <- full[apply(pvals, 1, function(r) sum(r < 5e-8) >=1),] # only keep variants sig in 1 or more metabolites; if want to be sig in at least 10 metabolites: sum(r < 5e-8) >=10)
pvals <- sig %>% select(matches("P_BOLT_LMM_"))
sig <- sig[apply(pvals, 1, function(r) sum(r < 1e-4) >=2),] # only keep variants sig in 1 or more metabolites; if want to be sig in at least 10 metabolites: sum(r < 5e-8) >=10)

mets <- gsub("BETA_","",colnames(sig %>% select(matches("BETA_"))))
sig <- distinct(sig)
sig_long <- data.table::melt(sig, measure.vars = patterns("^BETA_", "^SE_", "^P_BOLT_LMM"),variable.name=c("Metabolite"), value.name=c("BETA", "SE", "P_BOLT_LMM"))
sig_long$Metabolite <- rep(mets,each=nrow(sig_long)/length(mets))
sig_long <- filter(sig_long,P_BOLT_LMM < 1e-4)
sig_cart <- inner_join(sig_long,sig_long, by=c("SNP","Genes"), suffix=c(".1", ".2"))
sig_cart <- filter(sig_cart,Metabolite.1!=Metabolite.2)
sig_cart <- filter(sig_cart,P_BOLT_LMM.1 < 5e-8 | P_BOLT_LMM.2 < 5e-8)
sig_cart <- sig_cart %>% mutate(beta_opp=ifelse(sign(BETA.1)==sign(BETA.2),"No","Yes"))

# Identify discordant variants
sig_lfsr <- inner_join(sig_cart,gencor_lfsr,by=c("Metabolite.1"="p1","Metabolite.2"="p2"))
sig_lfsr <- sig_lfsr %>% mutate(contradict_gencor=ifelse((beta_opp=="Yes" & rg>0) | (beta_opp=="No" & rg<0),"Yes", "No"))

con <- left_join(filter(sig_lfsr,contradict_gencor=="No"),an,by=c("SNP")) %>% select(SNP,CHROM,BP,Assigned_Gene,Assigned_Gene_Type,matches(".1"),matches(".2$"),rg,se,z,p,lfsr)
confilt <- con %>% filter(Assigned_Gene_Type=="Pathway_Relevant_Enzyme") %>% arrange(desc(Assigned_Gene)) %>% select(SNP,Assigned_Gene,Assigned_Gene_Type,Metabolite.1,Metabolite.2)
confilt <- as.data.frame(confilt %>% mutate(Relative_Location="Not_Between",Temp=paste0(Assigned_Gene,Metabolite.1,Metabolite.2)))
confilt["Relative_Location"][confilt["Temp"]=="PCCBGlyVal"] <- "Between"
confilt["Relative_Location"][confilt["Temp"]=="C12orf5AlaGlucose"] <- "Between"
confilt <- confilt %>% select(-Temp)

v <- bind_rows(dfilt %>% mutate(Variant_Classification="Discordant"),confilt %>% mutate(Variant_Classification="Concordant"))

OUTPUT_FILE="discorvar_rellocation.tsv"
write.table(v, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

# for dis and con vars annotated with pathway relevant enzyme genes:
# num concordant with all met pairs "not between" = 12, num concordant with at least one met pair "between" = 2, num discordant with at least one met pair "not between" = 1, num discordant with all met pairs "between" = 5
fisher.test(matrix(c(12,2,1,5), nrow=2))
