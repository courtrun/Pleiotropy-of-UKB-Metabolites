# Courtney Smith - Nightingale Analysis - BOLT-LMM - Local gen cor analyses
# Goal of script: Table combining results across dif methods for local gen cor analysis of Alanine and Glutamine

library(dplyr)

# HE regression
he <- data.table::fread("all_hereg_regressions.tsv")
he <- filter(he,V1=="Ala"&V2=="Gln") %>% select(-V3,-V5,-V9)
colnames(he) <- c("Trait1","Trait2","Pathway","HE_estimate","HE_analytical_se","HE_jackknife_se")
he <- as.data.frame(he)
he["Pathway"][he["Pathway"]=="Coding"] <- "ALL_GENES"

# Fligner
f <- data.table::fread("all_fligner_variance_tests.tsv")
f <- filter(f,g1=="Ala"&g2=="Gln") %>% select(-parameter,-method) %>% rename(Trait1=g1,Trait2=g2,Pathway=pathway,Flinger_statistic=statistic,Flinger_P=p.value)
f <- as.data.frame(f)
f["Pathway"][f["Pathway"]=="Coding"] <- "ALL_GENES"

# HESS
h <- data.table::fread("pathway_group_tf_transporter_minmaf05_Ala_Gln_table.tsv")
h <- h %>% select(V7,V8,V2,V6,V17,V18,V19) %>% rename(Trait1=V7,Trait2=V8,Pathway=V2,HESS_LDblocks=V6,HESS_gencor=V17,HESS_se=V18,HESS_P=V19)
h <- as.data.frame(h)
h["Pathway"][h["Pathway"]=="GENOME"] <- "ALL_GENES"

# sLDSC
l <- data.table::fread("all_enrich_trivialmet_tabs.tsv")
colnames(l) <- c("Pathway","Trait1","Trait2","sLDSC_rg","sLDSC_se","sLDSC_z","sLDSC_p","sLDSC_h2_obs",
                "sLDSC_h2_obs_se","sLDSC_h2_int","sLDSC_h2_int_se","sLDSC_gcov_int","sLDSC_gcov_int_se")
l <- as.data.frame(l)
l["Pathway"][l["Pathway"]=="ANYGENE"] <- "ALL_GENES"

# BOLT-REML
b <- data.table::fread("bolt_reml_Ala_Gln_estimates.tsv") %>% mutate(Trait1="Ala",Trait2="Gln")
colnames(b) <- c("Pathway","REML_gencor","REML_se","Trait1","Trait2")
b <- as.data.frame(b)
b["Pathway"][b["Pathway"]=="Coding"] <- "ALL_GENES"

tb <- left_join(he,f,by=c("Trait1","Trait2","Pathway"))
tb <- left_join(tb,h,by=c("Trait1","Trait2","Pathway"))
tb <- left_join(tb,l,by=c("Trait1","Trait2","Pathway"))
tb <- left_join(tb,b,by=c("Trait1","Trait2","Pathway"))
keep_pathways <- c("ALL_GENES", "CUSTOM_TF","CUSTOM_TRANSPORTERS","CUSTOM_AA", "CUSTOM_KETONE", "CUSTOM_BCAA", "CUSTOM_GLYCOLYSIS_GLUCONEOGENESIS_TCA", "CUSTOM_UREAMET")
tb <- filter(tb,Pathway %in% keep_pathways)

FILE="combined_localgencor_methods.tsv"
write.table(tb,FILE, quote=F, sep="\t", row.names=F, col.names=T)
