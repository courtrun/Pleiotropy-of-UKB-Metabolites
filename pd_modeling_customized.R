# Courtney Smith - Nightingale Analysis - BOLT-LMM - Pathway Diagram Modeling
# Goal of script: Project a given variant's sumstats info onto pathway diagram - customized background

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

options(ggrepel.max.overlaps = Inf)

# Load in matrix
mGWAS <- data.table::fread("mGWAS_rigidmetabolites_sig_pruned.tsv")
cat <- data.table::fread("biomarkers_pathwayorder.csv") %>% select(Name_in_TSV,Color_Group)
lipids <- filter(cat,Color_Group=="Lipid_Particle")$Name_in_TSV
mGWAS <- select(mGWAS,!matches(lipids))
an <- data.table::fread("sigprunedmetabolites_filtfa_annotated_full.snplist.tsv")
an <- distinct(an %>% group_by(SNP) %>% arrange(Core_Gene) %>% select(-CHROM,-BP,-Pathway)) %>% summarise(Genes=paste(Core_Gene,collapse=","))
mGWAS_an <- left_join(mGWAS,an,by=c("SNP"))

# assign positions of respective metabolites
positions_dt <- data.frame(x=c(26,18,58.5,8,43,62,44.5,80,53,75,67,59.5,70,81,9.5,15), y= c(98.5,51,31.5,51,40.5,73,61.5,28.5,9,63.5,68.5,55.5,17.5,17.5,7.75,13.75), Metabolite = c("Glucose","Pyruvate","Gln","Lactate","Citrate","Gly","Ala","Acetoacetate","His","Val","Ile","Leu","Acetone","bOHbutyrate","Phe","Tyr"))

# Prepare matrix
mets <- gsub("BETA_","",colnames(mGWAS_an %>% select(matches("BETA_"))))
mGWAS_an_long <- data.table::melt(mGWAS_an, measure.vars = patterns("^BETA_", "^SE_", "^P_BOLT_LMM"),variable.name=c("Metabolite"), value.name=c("BETA", "SE", "P_BOLT_LMM"))
mGWAS_an_long$Metabolite <- rep(mets,each=nrow(mGWAS_an_long)/length(mets))

# combined positions with data to new data frame
mGWAS_an_positions <- left_join(mGWAS_an_long,positions_dt,by=c("Metabolite"))

pd_modeling <- function(var,my_gene_annotation){
# load in image to overlay data on
filename=paste0("rigid_pathways_",my_gene_annotation,".PNG")
r <- png::readPNG(filename)
rg <- grid::rasterGrob(r, width=unit(1,"npc"), height=unit(1,"npc"))

# filter to just var info
var_info <- filter(mGWAS_an_positions,SNP==var) %>% mutate(Gene_annotation=my_gene_annotation)
if (my_gene_annotation=="PCCB" | my_gene_annotation=="GCKR"){
l <- data.table::fread("INT_LDL_direct_adjstatins_all.glm.linear.filtered.maf001.info03.tsv") # ukbiobank
h <- data.table::fread("INT_HDL_cholesterol_all.glm.linear.filtered.maf001.info03.tsv") # ukbiobank
f <- data.table::fread("Total_FA/european_phe_and_remove.biomarkers.imp_stats.filt_sigclumped.gz")
tg <- data.table::fread("Total_TG/european_phe_and_remove.biomarkers.imp_stats.filt_sigclumped.gz")
l <- filter(l,SNP==var) %>% mutate(Metabolite="LDL_C",Genes.x="",Genes.y="",x=65,y=45,Gene_annotation=my_gene_annotation,P_BOLT_LMM=as.numeric(P),BETA=-BETA)  %>% # flip sign to align alleles
      select(SNP,Genes.x,Genes.y,Metabolite,BETA,SE,P_BOLT_LMM,x,y,Gene_annotation)
h <- filter(h,SNP==var) %>% mutate(Metabolite="HDL_C",Genes.x="",Genes.y="",x=69,y=43,Gene_annotation=my_gene_annotation,P_BOLT_LMM=as.numeric(P),BETA=-BETA)  %>% # flip sign to align alleles
      select(SNP,Genes.x,Genes.y,Metabolite,BETA,SE,P_BOLT_LMM,x,y,Gene_annotation)
f <- filter(f,SNP==var) %>% mutate(Metabolite="Total_FA",Genes.x="",Genes.y="",x=61.5,y=54.5,Gene_annotation=my_gene_annotation,P_BOLT_LMM=as.numeric(P_BOLT_LMM)) %>% select(SNP,Genes.x,Genes.y,Metabolite,BETA,SE,P_BOLT_LMM,x,y,Gene_annotation)
tg <- filter(tg,SNP==var) %>% mutate(Metabolite="Total_TG",Genes.x="",Genes.y="",x=45.5,y=35.5,Gene_annotation=my_gene_annotation,P_BOLT_LMM=as.numeric(P_BOLT_LMM)) %>% select(SNP,Genes.x,Genes.y,Metabolite,BETA,SE,P_BOLT_LMM,x,y,Gene_annotation)
var_info <- bind_rows(var_info,l,h,f,tg)}

min_lim=-0.07
max_lim=-min_lim
if (my_gene_annotation=="CPS1"){
var_info <- as.data.frame(var_info)
var_info["BETA"][var_info["Metabolite"]=="Gly"] <- -0.07
}
# PLOTS for each variant, colored by beta effect size and shape indicating significance
ggplot(var_info,aes(x,y,fill=BETA)) + # ,alpha=P_BOLT_LMM<1e-4
  annotation_custom(rg) + geom_point(size=3,color=ifelse(var_info$P_BOLT_LMM<1e-4,"black","gray"),pch=21) + #scale_shape_manual(values=c(22,21))
  scale_x_continuous(expand=c(0,0), lim=c(0,100)) +
  scale_y_continuous(expand=c(0,0), lim=c(0,100)) +
  scale_fill_gradientn(colors=c("darkblue","blue","white", "orange","darkorange"),limits=c(min_lim,max_lim))+ # ,limits=c(-max_beta, max_beta) # colors=c("purple","blue","cyan","green","white", "yellow","orange","red","darkred")
  #geom_label_repel(aes(label=paste0(Metabolite)), size=4, fill="white", min.segment.length=unit(0,'lines'),force=10,segment.size=0.25, show.legend = FALSE)+
  theme_void() +
  theme(aspect.ratio = nrow(r)/ncol(r)) +labs(fill="Beta",title=paste0("Effect of ",var," (",my_gene_annotation,") on Relevant Metabolites"))
print(ggsave(paste0("pd_modeling_",var,".pdf")))
}

an2 <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation.tsv")


genes <- c("PCCB")
for (gene in genes) {
var=filter(an2,Gene==gene)$SNP
if (length(var)>1) {
for (v in var) {
var=v
my_gene_annotation=gene
pd_modeling(var,my_gene_annotation)
}}
my_gene_annotation=gene
pd_modeling(var,my_gene_annotation)
}

var="rs77010315"
my_gene_annotation="SLC36A2"
pd_modeling(var,my_gene_annotation)

var="rs370014171"
my_gene_annotation="PDPR"
pd_modeling(var,my_gene_annotation)
