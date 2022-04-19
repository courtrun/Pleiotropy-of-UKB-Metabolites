# Courtney Smith - Nightingale Analysis - BOLT-LMM - Discordant variant analysis
# Goal of script: Bar plot of gene type by variant type enrichment

thres1 = 1e-4
thres2 = 1e-4

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
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation.tsv")

# filter to pairs sig dif than 0
filt <- ashr::ash(gc$rg,gc$se)
lfsr <- as.data.frame(filt$result) %>% select(lfsr)
gencor_lfsr <- cbind(gc,lfsr)
gencor_lfsr <- filter(gencor_lfsr,lfsr<0.005) # 0.005 because that's less than 5 in 1000 and theres's about 120 correlations

# Identify variants sig in at least one rigid metabolite 5e-8 and 1e-4 in another
pvals <- full %>% select(matches("P_BOLT_LMM_"))
sig <- full[apply(pvals, 1, function(r) sum(r < thres1) >=1),] # only keep variants sig in 1 or more metabolites; if want to be sig in at least 10 metabolites: sum(r < 5e-8) >=10)
pvals <- sig %>% select(matches("P_BOLT_LMM_"))
sig <- sig[apply(pvals, 1, function(r) sum(r < thres2) >=2),] # only keep variants sig in 1 or more metabolites; if want to be sig in at least 10 metabolites: sum(r < 5e-8) >=10)

mets <- gsub("BETA_","",colnames(sig %>% select(matches("BETA_"))))
sig <- distinct(sig)
sig_long <- data.table::melt(sig, measure.vars = patterns("^BETA_", "^SE_", "^P_BOLT_LMM"),variable.name=c("Metabolite"), value.name=c("BETA", "SE", "P_BOLT_LMM"))
sig_long$Metabolite <- rep(mets,each=nrow(sig_long)/length(mets))
sig_long <- filter(sig_long,P_BOLT_LMM < thres2)
sig_cart <- inner_join(sig_long,sig_long, by=c("SNP","Genes"), suffix=c(".1", ".2"))
sig_cart <- filter(sig_cart,Metabolite.1!=Metabolite.2)
sig_cart <- filter(sig_cart,P_BOLT_LMM.1 < thres1 | P_BOLT_LMM.2 < thres1)
sig_cart <- sig_cart %>% mutate(beta_opp=ifelse(sign(BETA.1)==sign(BETA.2),"No","Yes"))

# Identify discordant variants
sig_lfsr <- inner_join(sig_cart,gencor_lfsr,by=c("Metabolite.1"="p1","Metabolite.2"="p2"))
sig_lfsr <- sig_lfsr %>% mutate(contradict_gencor=ifelse((beta_opp=="Yes" & rg>0) | (beta_opp=="No" & rg<0),"Yes", "No"))

# How many variants are eligible (bc they are sig in at least one rigid metabolite 5e-8 and 1e-4 in another for pairs of mets that have gen cor sig dif than 0)
elig <- unique(sig_lfsr$SNP)

# What is the breakdown by gene type for discordant variants in at least one pair (note, some of these variants are annotated with the same gene...)
dis <- unique(filter(sig_lfsr,contradict_gencor=="Yes")$SNP)
con <- unique(filter(sig_lfsr,contradict_gencor=="No")$SNP)

### Run analysis
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation.tsv")
an$Gene_Type <- ifelse(an$Gene_Type == "Pathway_Relevant_Enzyme" | an$Gene_Type == "Transporter", "Enzyme or\nTransporter", an$Gene_Type)
an$Gene_Type <- ifelse(an$Gene_Type != "Enzyme or\nTransporter", "Not Enzyme\nor Transporter", an$Gene_Type)
an <- left_join(an,an %>% group_by(Gene_Type) %>% count() %>% rename(num_in_genetype=n),by=c("Gene_Type"))
an <- an %>% mutate(Variant_Type=ifelse(SNP %in% dis,"Discordant",ifelse(SNP %in% con,"Concordant","Neither")))
an <- left_join(an,an %>% group_by(Variant_Type) %>% count() %>% rename(num_in_vartype=n),by=c("Variant_Type"))

# combine transporter and enzyme
ansum <- left_join(an %>% group_by(Gene_Type,Variant_Type) %>% count(),an %>% group_by(Gene_Type,Variant_Type) %>% summarize(num_in_vartype=mean(num_in_vartype), num_in_genetype=mean(num_in_genetype)),by=c("Gene_Type","Variant_Type"))
ansum <- ansum %>% mutate(n_per_vt=n/num_in_vartype,n_per_gt=n/num_in_genetype,se_vt = sqrt(n*(num_in_vartype-n)/(num_in_vartype^3)),se_gt = sqrt(n*(num_in_genetype-n)/(num_in_genetype^3)))# se_vt=sqrt(n)/num_in_vartype,se_gt=sqrt(n)/num_in_genetype)

library("ggplot2")

# Fraction of variants of a given variant type in each gene type / # variants in that variant type vs gene type, fill by var type
ggplot(filter(ansum,Variant_Type!="Neither"),aes(y=Variant_Type,x=n_per_vt,fill=Gene_Type))+ geom_bar(stat="identity",  position=position_dodge2(reverse = TRUE),width=0.6)+
	labs(y="Variant Type",x="Fraction of Variants",fill="Gene Type")+
	scale_fill_manual(values=c("#6CFA9C","black"))+ # #6CFA9C, #85e085
	geom_errorbar(aes(xmin=n_per_vt-se_vt*2, xmax=n_per_vt+se_vt*2,y=Variant_Type), width=.6,position=position_dodge2(width = .9,reverse = TRUE, padding = .6)) +
	theme_classic()+theme(text=element_text(size=16)) # axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
ggsave("barplot_numbervarpervartype_twocat_disvscon.pdf",height=4,width=9)

d_et <- filter(ansum,Gene_Type=="Enzyme or\nTransporter"&Variant_Type=="Discordant")$n # num discordant enzymes/transporters
d_net <- filter(ansum,Gene_Type=="Not Enzyme\nor Transporter"&Variant_Type=="Discordant")$n # num discordant not enzymes/transporters
d_total <- filter(ansum,Gene_Type=="Enzyme or\nTransporter"&Variant_Type=="Discordant")$num_in_vartype # total num discordant
c_et <- filter(ansum,Gene_Type=="Enzyme or\nTransporter"&Variant_Type=="Concordant")$n # num Concordant enzymes/transporters
c_net <- filter(ansum,Gene_Type=="Not Enzyme\nor Transporter"&Variant_Type=="Concordant")$n # num Concordant not enzymes/transporters
c_total <- filter(ansum,Gene_Type=="Enzyme or\nTransporter"&Variant_Type=="Concordant")$num_in_vartype # total num Concordant
total_et <- filter(ansum,Gene_Type=="Enzyme or\nTransporter"&Variant_Type=="Discordant")$num_in_genetype # total num enzymes/transporters
total_net <- filter(ansum,Gene_Type=="Not Enzyme\nor Transporter"&Variant_Type=="Discordant")$num_in_genetype # total num not enzymes/transporters

# 95% CI calculated by binomial sampling variance
binom.test(d_et,d_total,total_et/(total_et+total_net)) # 0.04 # 116/(116+97) that's the ratio of enzymes/transporters vs not so gives how much would expect under null
binom.test(c_et,c_total,total_et/(total_et+total_net)) # 0.48

d <- fisher.test(matrix(c(d_et,d_net,total_et-d_et,total_net-d_net),nrow=2)) # odds ratio 4.09 and p value 0.034
c <- fisher.test(matrix(c(c_et,c_net,total_et-c_et,total_net-c_net),nrow=2)) # odds ratio 0.75, p-value 0.42
