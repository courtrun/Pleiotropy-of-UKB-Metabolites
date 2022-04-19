# Courtney Smith - Nightingale Analysis - BOLT-LMM - Metabolite Group Assignments
# Goal of script: Assign each Assigned_Gene that is a met hit to a (or multiple) metabolite group

library(dplyr)

an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")
full <- data.table::fread("mGWAS_rigidmetabolites_sig_pruned.tsv") %>% select(-Genes)
cat <- data.table::fread("biomarkers_pathwayorder.csv") %>% select(Name_in_TSV,Color_Group)
sig_thres = 1e-4 # 1e-4, 5e-8

# Add BCAA category
cat <- as.data.frame(cat)
BCAA_list <- c("Leu","Ile","Val")
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[1]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[2]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[3]] <- "BCAA"

# For each trait, count how many sig variant associations that enzyme has with that trait
mets <- gsub("BETA_","",colnames(full %>% select(matches("BETA_"))))
mets_long <- data.table::melt(full, measure.vars = patterns("^BETA_", "^SE_", "^P_BOLT_LMM"),variable.name=c("Metabolite"), value.name=c("BETA", "SE", "P_BOLT_LMM"))
mets_long$Metabolite <- rep(mets,each=nrow(mets_long)/length(mets))

# Add Assigned_Gene annotations and trait group assignments
mets_long <- left_join(mets_long,an,by=c("SNP"))
mets_long <- left_join(mets_long,cat,by=c("Metabolite"="Name_in_TSV"))

# Filter to sig variant-trait associations for this Assigned_Gene type
e <- mets_long %>% filter(Assigned_Gene_Type == "Pathway_Relevant_Enzyme") #
e <- filter(e,P_BOLT_LMM<=sig_thres)
ge <- as.data.frame(as.data.frame(e %>% group_by(Assigned_Gene,Color_Group) %>% count()) %>% group_by(Assigned_Gene) %>% top_n(1, n))
ge_all <- as.data.frame(as.data.frame(e %>% group_by(Assigned_Gene,Color_Group) %>% count()))

# Filter to sig variant-trait associations for this Assigned_Gene type
t <- mets_long %>% filter(Assigned_Gene_Type == "Transporter") #
t <- filter(t,P_BOLT_LMM<=sig_thres)
gt <- as.data.frame(as.data.frame(t %>% group_by(Assigned_Gene,Color_Group) %>% count()) %>% group_by(Assigned_Gene) %>% top_n(1, n))
filter(gt,Color_Group=="Ketone")
filter(gt,Color_Group=="Glycolysis")
filter(gt,Color_Group=="Amino_Acid")
filter(gt,Color_Group=="BCAA")
gt_all <- as.data.frame(as.data.frame(t %>% group_by(Assigned_Gene,Color_Group) %>% count()))

# Filter to sig variant-trait associations for this Assigned_Gene type
tf <- mets_long %>% filter(Assigned_Gene_Type == "Transcription") #
tf <- filter(tf,P_BOLT_LMM<=sig_thres)
gtf <- as.data.frame(as.data.frame(tf %>% group_by(Assigned_Gene,Color_Group) %>% count()) %>% group_by(Assigned_Gene) %>% top_n(1, n))
filter(gtf,Color_Group=="Ketone")
filter(gtf,Color_Group=="Glycolysis")
filter(gtf,Color_Group=="Amino_Acid")
filter(gtf,Color_Group=="BCAA")
gtf_all <- as.data.frame(as.data.frame(tf %>% group_by(Assigned_Gene,Color_Group) %>% count()))

# Save as table combined across all - one assignment per gene
ge <- ge %>% mutate(Assigned_Gene_Type="Pathway_Relevant_Enzyme") %>% select(Assigned_Gene,Assigned_Gene_Type,Color_Group,n) %>% rename(Assigned_Metabolite_Group=Color_Group,Metabolite_Group_Associations=n)
gt <- gt %>% mutate(Assigned_Gene_Type="Transporter")  %>% select(Assigned_Gene,Assigned_Gene_Type,Color_Group,n)%>% rename(Assigned_Metabolite_Group=Color_Group,Metabolite_Group_Associations=n)
gtf <- gtf %>% mutate(Assigned_Gene_Type="Transcription") %>% select(Assigned_Gene,Assigned_Gene_Type,Color_Group,n) %>% rename(Assigned_Metabolite_Group=Color_Group,Metabolite_Group_Associations=n)
g <- bind_rows(ge,gt,gtf)
g$Assigned_Metabolite_Group <- ifelse(g$Assigned_Metabolite_Group=="Amino_Acid","Other Amino Acids",ifelse(g$Assigned_Metabolite_Group=="Ketone","Ketone Bodies",g$Assigned_Metabolite_Group))

OUTPUT_FILE="gene_metgroup_assignments.tsv"
write.table(g, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

# Save as table combined across all - all assignments per gene
ge_all <- ge_all %>% mutate(Assigned_Gene_Type="Pathway_Relevant_Enzyme") %>% select(Assigned_Gene,Assigned_Gene_Type,Color_Group,n) %>% rename(Assigned_Metabolite_Group=Color_Group,Metabolite_Group_Associations=n)
gt_all <- gt_all %>% mutate(Assigned_Gene_Type="Transporter")  %>% select(Assigned_Gene,Assigned_Gene_Type,Color_Group,n)%>% rename(Assigned_Metabolite_Group=Color_Group,Metabolite_Group_Associations=n)
gtf_all <- gtf_all %>% mutate(Assigned_Gene_Type="Transcription") %>% select(Assigned_Gene,Assigned_Gene_Type,Color_Group,n) %>% rename(Assigned_Metabolite_Group=Color_Group,Metabolite_Group_Associations=n)
g_all <- bind_rows(ge_all,gt_all,gtf_all)
g_all$Assigned_Metabolite_Group <- ifelse(g_all$Assigned_Metabolite_Group=="Amino_Acid","Other Amino Acids",ifelse(g_all$Assigned_Metabolite_Group=="Ketone","Ketone Bodies",g_all$Assigned_Metabolite_Group))

OUTPUT_FILE="gene_metgroup_assignments_all.tsv"
write.table(g_all, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)
