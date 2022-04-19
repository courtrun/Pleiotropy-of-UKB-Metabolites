# Courtney Smith - Nightingale Analysis - BOLT-LMM - Local rho
# Goal of script: Make pathway regions for local rho HESS analysis

### Prepare combined_pathway.txt file for automated kegg rigid pathways
library(dplyr)

c1 <- data.table::fread("REACTOME_REGULATION_OF_PYRUVATE_DEHYDROGENASE_PDH_COMPLEX.bed")
c2 <- data.table::fread("GO_REGULATION_OF_HEXOKINASE_ACTIVITY.bed")
c3 <- data.table::fread("GO_KETONE_BODY_METABOLIC_PROCESS.bed")
c4 <- data.table::fread("KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION.bed")
c5 <- data.table::fread("KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM.bed")
c6 <- data.table::fread("KEGG_HISTIDINE_METABOLISM.bed")
c7 <- data.table::fread("KEGG_TYROSINE_METABOLISM.bed")
c8 <- data.table::fread("KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM.bed")
c9 <- data.table::fread("KEGG_PYRUVATE_METABOLISM.bed")
c10 <- data.table::fread("KEGG_CITRATE_CYCLE_TCA_CYCLE.bed")
c11 <- data.table::fread("KEGG_GLYCOLYSIS_GLUCONEOGENESIS.bed")

c_all <- bind_rows(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)
c_all <- c_all %>% mutate(V2=V2+100000,V3=V3-100000) # adjust back to just gene boundaries
c_all$V4 <- gsub("_([^_]*)$", c_all$V4, replacement=" \\1")

p <- tidyr::separate(c_all,"V4",into=c("Pathway","Gene"),sep=" ") %>% arrange(V1,V2)

OUTPUT_FILE="combined_pathway.txt"
write.table(p, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

### Prepare custom rigid file
library(dplyr)

g <- data.table::fread("allpathways_msig_within100kbOfGenes_wpathways.bed")
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")

# Make gene list
genes <-c("ACO1","PYCR2","PYCRL","PYCR1","ALDH4A1","ALDH18A1","FBP1","FBP2","PCK2",
"SUCLG1","SUCLG2","SDHA","SDHB","SDHC","SDHD","SDS","AMDHD1","FTCD","PSAT1","ACO2",
"AUH","MCCC1","HMGCL","ASS1","ALDOA","PFKM","GPT","H6PD","PHGDH","BCAT2","GCKR","PFKFB3",
"SHMT1","GPI","G6PD","PKM2","ENO2","PGAM1","ENO1","FH","CS","GCSH","HMGCS1","HMGCS2","GCK",
"PKLR","LDHA","LDHB","LDHC","PDHX","PDK1","PDK2","IDH1","MTHFS","ARG1","CKMT1A","MSRA","GCH1",
"GAPDH","BDH2","ARG2","HOGA1","ASL","PC","IVD","GLUL","ACAT1","ACAT2","OTC","SHMT2","HPD",
"PCK1","ACAA1","SUCLA2","PGK1","PGK2","MDH1","MDH2","GATM","TDO2","IDO1","IDO2","AFMID",
"KMO","KYNU","MTHFD1","AMT","OAT")
pre <- filter(an,Assigned_Gene_Type=="Pathway_Relevant_Enzyme")$Assigned_Gene
bcaa <- data.table::fread("KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION.bed")
bcaa$V4 <- gsub("_([^_]*)$", bcaa$V4, replacement=" \\1")
bcaa <- tidyr::separate(bcaa,"V4",into=c("Pathway","V4"),sep=" ") %>% arrange(V1,V2)
skip <- c("OXCT2","HMGCL","HMGCS1","HMGCS2","ACAT1","ACAT2","OXCT1","DLD","AGXT2","ABAT","AACS")
bcaa <- filter(bcaa,!(V4 %in% skip))$V4
glist <- unique(c(genes,pre,bcaa))
skip <- c("ACOT11","FDFT1","CKMT1A","HOGA1","H6PD","G6PD","FADS1")
glist <- setdiff(glist,skip)

# Make file with just those genes
r <- filter(g,V4 %in% glist)
c <- distinct(r %>% select(-V5)) %>% mutate(V5="CUSTOM_RIGID") %>% arrange(V1,V2) %>% select(V1,V2,V3,V5,V4)
c <- c %>% mutate(V2=V2+100000,V3=V3-100000)

# Assign individual pathways
ketone <- c("ACAT1","ACAT2","OXCT1", "BDH2", "HMGCS1", "HMGCS2", "HMGCL")
urea.met <- c("ALDH1L1","MTHFS", "MTHFD1", "AMT", "MTHFR", "MSRA", "GCSH", "GLDC", "GATM", "CKMT1A", "CPS1", "OTC", "ARG1", "ARG2", "ASS1", "ASL")
gly <- c("GCK","GCKR","GPI","PFKFB2","PFKFB3","C12orf5","PFKL","PFKM","PFKP","ALDOA","GAPDH","PGK1","PGK2","PGM2L1","PGAM1","ENO1","ENO2","PKLR","LDHA","LDHB","LDHC","PDHX","PDPR","PDK1","PDK2","PDK4","PCK1","PCK2","PC","FBP1","FBP2","G6PC","G6PC2","ACAA1","CS","ACO1","ACO2","IDH1","DLST","SUCLA2","SUCLG1","SUCLG2","SDHA","SDHB","SDHC","SDHD","FH","MDH1","MDH2")
bcaa <- c(bcaa,"ECHDC1")
aa <- unique(filter(c,!(V4 %in% c(ketone,urea.met,gly,bcaa)))$V4)

c <- c %>% mutate(V5=ifelse(V4 %in% aa,"CUSTOM_AA",ifelse(V4 %in% gly,"CUSTOM_GLYCOLYSIS_GLUCONEOGENESIS_TCA",ifelse(V4 %in% ketone,"CUSTOM_KETONE",ifelse(V4 %in% urea.met,"CUSTOM_UREAMET",ifelse(V4 %in% bcaa,"CUSTOM_BCAA","NA"))))))
c <- filter(c,!(V4=="SDHD"&V3=="112064528"))

OUTPUT_FILE="CUSTOM_RIGID.bed"
write.table(c, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=F)

### Make pathway list with all genes and their pathway (no extension from normal gene boundaries; if was extended then I unextended)
library(dplyr)

c1 <- data.table::fread("combined_pathway.txt")
c2 <- data.table::fread("CUSTOM_RIGID.bed")
c <- bind_rows(c1,c2)

c4 <- data.table::fread("REACTOME_CHOLESTEROL_BIOSYNTHESIS.bed")
c5 <- data.table::fread("REACTOME_CHYLOMICRON_MEDIATED_LIPID_TRANSPORT.bed")
c6 <- data.table::fread("KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS.bed")
c7 <- data.table::fread("GO_NITROGEN_CYCLE_METABOLIC_PROCESS.bed")
c8 <- data.table::fread("REACTOME_METABOLISM_OF_POLYAMINES.bed")
c9 <- data.table::fread("GO_METHIONINE_METABOLIC_PROCESS.bed")

c_all <- bind_rows(c4,c5,c6,c7,c8,c9)
c_all <- c_all %>% mutate(V2=V2+100000,V3=V3-100000)
c_all$V4 <- gsub("_([^_]*)$", c_all$V4, replacement=" \\1")

c_all <- tidyr::separate(c_all,"V4",into=c("V4","V5"),sep=" ")
p <- bind_rows(c,c_all)
colnames(p) <- c("CHROM","POS_start","POS_end","Pathway","Gene")

AMINO_ACID_ALL <- c("KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM","KEGG_HISTIDINE_METABOLISM","KEGG_TYROSINE_METABOLISM","KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM")
AMINO_ACID_noBCAA <- c("KEGG_HISTIDINE_METABOLISM","KEGG_TYROSINE_METABOLISM","KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM")
BCAA <- c("KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM")
GLYCOLYSIS_GLUCONEOGENESIS_TCA <- c("KEGG_PYRUVATE_METABOLISM","KEGG_CITRATE_CYCLE_TCA_CYCLE","KEGG_GLYCOLYSIS_GLUCONEOGENESIS","REACTOME_REGULATION_OF_PYRUVATE_DEHYDROGENASE_PDH_COMPLEX","GO_REGULATION_OF_HEXOKINASE_ACTIVITY")
KETONE <- c("GO_KETONE_BODY_METABOLIC_PROCESS")
LIPID <- c("REACTOME_CHOLESTEROL_BIOSYNTHESIS","REACTOME_CHYLOMICRON_MEDIATED_LIPID_TRANSPORT","KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS")
UREA_METHIONINE <- c("GO_NITROGEN_CYCLE_METABOLIC_PROCESS","GO_METHIONINE_METABOLIC_PROCESS","REACTOME_METABOLISM_OF_POLYAMINES")

p <- p %>% mutate(Pathway_Group=ifelse(Pathway %in% AMINO_ACID_ALL,"AMINO_ACID_ALL",ifelse(Pathway %in% GLYCOLYSIS_GLUCONEOGENESIS_TCA,"GLYCOLYSIS_GLUCONEOGENESIS_TCA",ifelse(Pathway %in% KETONE,"KETONE",ifelse(Pathway %in% LIPID,"LIPID",ifelse(Pathway %in% UREA_METHIONINE,"UREA_METHIONINE",ifelse(Pathway %in% AMINO_ACID_noBCAA,"AMINO_ACID_noBCAA",ifelse(Pathway %in% BCAA,"BCAA",ifelse(startsWith(Pathway,"CUSTOM"),"CUSTOM_RIGID","NA")))))))))

OUTPUT_FILE="combined_pathway_wlipids_ureamet.txt"
write.table(p, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

### Add transporters and TFs to combined_pathway_wlipids_ureamet.txt (no extension from normal gene boundaries; if was extended then I unextended)
library(dplyr)

c <- data.table::fread("combined_pathway_wlipids_ureamet.txt")

# Make custom transporter and TF pathways
g <- data.table::fread("allpathways_msig_within100kbOfGenes_wpathways.bed")
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")
tf <- unique(filter(an,Assigned_Gene_Type=="Transcription")$Assigned_Gene)
c1 <- distinct(filter(g,V4 %in% tf) %>% select(V1,V2,V3,V4)) %>% mutate(Pathway="CUSTOM_TF") %>% select(V1,V2,V3,Pathway,V4)
tr <- unique(filter(an,Assigned_Gene_Type=="Transporter")$Assigned_Gene)
c2 <- distinct(filter(g,V4 %in% tr) %>% select(V1,V2,V3,V4)) %>% mutate(Pathway="CUSTOM_TRANSPORTERS") %>% select(V1,V2,V3,Pathway,V4)
c_12 <- bind_rows(c1,c2)
c_12 <- c_12 %>% mutate(V2=V2+100000,V3=V3-100000)
setdiff(tf,c1$V4) # ZNF764 chr16:30,565,085-30,569,819
setdiff(tr,c2$V4) # MPC1 chr6:166,778,407-166,796,486
extra_genes <- data.frame(V1=c("chr16","chr6"),V2=c(30565085,166778407),V3=c(30569819,166796486),Pathway=c("CUSTOM_TF","CUSTOM_TRANSPORTERS"),V4=c("ZNF764","MPC1"))
c_12 <- bind_rows(c_12,extra_genes)

# Add premade transporter and TF pathways
c3 <- data.table::fread("GO_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER.bed")
c4 <- data.table::fread("REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT.bed")

c_all <- bind_rows(c3,c4)
c_all <- c_all %>% mutate(V2=V2+100000,V3=V3-100000)
c_all$V4 <- gsub("_([^_]*)$", c_all$V4, replacement=" \\1")
c_all <- tidyr::separate(c_all,"V4",into=c("Pathway","V4"),sep=" ")

p <- bind_rows(c_12,c_all) %>% mutate(Pathway_Group=Pathway)
colnames(p) <- c("CHROM","POS_start","POS_end","Pathway","Gene","Pathway_Group")

a <- bind_rows(c,p)

OUTPUT_FILE="combined_pathway_wlipids_ureamet_tfandtransporters.txt"
write.table(a, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)
