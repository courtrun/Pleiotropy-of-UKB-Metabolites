# Courtney Smith - Nightingale Analysis - BOLT-LMM - Manhattan plots
# Goal of script: Plot manhattan plot for metabolites by biochemial group

# Load the library
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(hash)

# Read in file with metabolite categories for coloring
INPUT_FILE_CAT="biomarkers_pathwayorder.csv"
cat <- data.table::fread(INPUT_FILE_CAT)
# Add BCAA category
cat <- as.data.frame(cat)
BCAA_list <- c("Leu","Ile","Val")
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[1]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[2]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[3]] <- "BCAA"

an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")

options(ggrepel.max.overlaps = Inf)

# Set color dict by group key
color_dict <- hash()
color_dict[["Ketone"]] <- "purple"
color_dict[["Lipid_Particle"]] <- "blue"
color_dict[["Amino_Acid"]] <- "orange"
color_dict[["Glycolysis"]] <- "red"
color_dict[["Fatty_Acid"]] <- "green"
color_dict[["Other"]] <- "pink"

# Set color dict by color key
color_dict[["purple"]] <- c("darkorchid4","darkorchid1")
color_dict[["blue"]] <- c("dodgerblue4","cornflowerblue")
color_dict[["orange"]] <- c("darkorange2","orange")
color_dict[["red"]] <- c("darkred","firebrick1")
color_dict[["green"]] <- c("forestgreen","chartreuse3")
color_dict[["pink"]] <- c("deeppink4","deeppink")

# Store function that plots and saves manhattan plot given any metabolite name
manh_plot <- function(mets) {

met_list <- list()

# Loop through and get sumstats for each metabolite
for (met in mets) {
FILE=paste0("metabolomeFirstPass.",met,".annotated.bed")
if (file.exists(FILE)){
# Read in file with annotated data and prepare it
hits <- data.table::fread(FILE)} else { next}
colnames(hits) <- c('CHROM','POS','POS_end','ID','BETA','SE','P_BOLT_LMM','CHR_gene','chromStart','chromEnd','GENE')
hits$CHROM <- gsub("chr", "", hits$CHROM) # remove "chr" prefix
hits$CHROM <- as.numeric(hits$CHROM)
hits <- hits %>% group_by(CHROM,POS,ID,BETA,SE,P_BOLT_LMM) %>% arrange(GENE) %>% summarize(GENES=paste(GENE,collapse=",")) # one row per variant, list genes
hits$P_BOLT_LMM <- as.numeric(hits$P_BOLT_LMM)
hits$P_BOLT_LMM <- ifelse(hits$P_BOLT_LMM==0,5e-324,hits$P_BOLT_LMM) # set any pvalues = 0 to 5e-324
hits <- hits %>% filter(-log10(P_BOLT_LMM)>3) # to make plotting faster, filter out the very low significance variants
# Add column with cummulative position in genome of each variant, find lowest P SNP and value
nCHR <- length(unique(hits$CHROM))
hits$BPcum <- NA
s <- 0
nbp <- c(0)
for (i in unique(hits$CHROM)){
nbp[i+1] <- max(hits[hits$CHROM == i,]$POS) + s # store max POS for each chrom and add to previous chrom length
s <- nbp[i+1]  # Calculate new cumulative position so far
}
hits <- hits %>% mutate(BPcum=POS+nbp[CHROM]) # POS in chrom + cummulative base pairs before that chrom
met_list[[met]] <- data.frame(hits,Metabolite=met)
}

group <- filter(cat,Name_in_TSV==met)$Color_Group
mets_hits <- bind_rows(met_list)

# Calculate axes locations
axis.set <- mets_hits %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(mets_hits$P_BOLT_LMM)))) + 2
sig <- 5e-8

mets_hits <- left_join(mets_hits,an,by=c("ID"="SNP"),suffix=c("",".y"))
mets_hits$Gene_label <- ifelse(is.na(mets_hits$Assigned_Gene),mets_hits$GENES,mets_hits$Assigned_Gene)
mets_hits$Gene_Type_full <- ifelse(is.na(mets_hits$Assigned_Gene),"Not Manually Annotated",mets_hits$Assigned_Gene_Type)

snp_keep_gene <- (as.data.table(mets_hits)[ , .SD[which.min(P_BOLT_LMM)], by = Assigned_Gene])$ID  # filter to lowest p for each gene
snp_keep <- (as.data.table(filter(mets_hits,ID %in% snp_keep_gene))[ , .SD[which.min(P_BOLT_LMM)], by = CHROM])$ID # then keep lowest P for each CHROM

mets_hits$is_annotate <- ifelse(mets_hits$ID %in% snp_keep & mets_hits$P_BOLT_LMM < 5e-8,"yes","no")

# if there are multiple variants with the same rsID to be annotated, only do the one with the lowest
mets_hits <- mets_hits %>% mutate(uniqID=paste0(ID,Metabolite))
snp_keep_id <- (as.data.table(filter(mets_hits,is_annotate=="yes"))[ , .SD[which.min(P_BOLT_LMM)], by = ID])$uniqID # first filter to just lowest P for each
mets_hits$is_annotate <- ifelse(mets_hits$uniqID %in% snp_keep_id,"yes","no")

# manually update gene label for a var
mets_hits <- as.data.frame(mets_hits)
mets_hits["Gene_label"][mets_hits["Gene_label"]=="MED31,SLC13A5"] <- "SLC13A5"

# if a variant was not a rigid met hit but is annotated with a rigid met gene, add gene type
mets_hits <- left_join(mets_hits,distinct(an %>% select(Assigned_Gene,Assigned_Gene_Type) %>% filter(Assigned_Gene!=".")) %>% rename(Gene_Type_rm=Assigned_Gene_Type),by=c("Gene_label"="Assigned_Gene"))
mets_hits$Gene_Type_full <- ifelse(mets_hits$Gene_Type_full=="Not Manually Annotated" & !is.na(mets_hits$Gene_Type_rm),mets_hits$Gene_Type_rm,mets_hits$Gene_Type_full)

dark_colors <- RColorBrewer::brewer.pal(n=8,name="Dark2")

# only color if is in line with labeled variant
mets_hits <- mets_hits %>% mutate(Gene_Type_full_color=ifelse(is_annotate=="yes",Gene_Type_full,"Not Manually Annotated"))

if (length(unique(mets_hits$Gene_Type_full_color))==7) {
color_list = c(dark_colors[2],dark_colors[7],"black",dark_colors[4],dark_colors[3],dark_colors[1],dark_colors[6])
} else if (length(unique(mets_hits$Gene_Type_full_color))==6) {
color_list = c(dark_colors[2],dark_colors[7],"black",dark_colors[4],dark_colors[1],dark_colors[6])
} else if (length(unique(mets_hits$Gene_Type_full_color))==5){
color_list = c(dark_colors[2],"black",dark_colors[4],dark_colors[1],dark_colors[6])
} else if (length(unique(mets_hits$Gene_Type_full_color))==4){
color_list = c(dark_colors[2],"black",dark_colors[4],dark_colors[1])
}

mets_hits["Assigned_Gene_Type"][mets_hits["Assigned_Gene_Type"]=="Transcription"] <- "TF"

# Plot manhattan plot
ggplot(mets_hits,aes(x=BPcum,y=-log10(P_BOLT_LMM),color=as.factor(Gene_Type_full_color)))+geom_point(alpha=0.75,size=0.5)+ #,color=as.factor(CHROM) # +geom_point(color=as.factor(hits$CHROM),alpha=0.75,size=0.5)+
geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + # add significance line
scale_x_continuous(label = axis.set$'CHROM', breaks = axis.set$center) + # label the chromosome number in the center of each chromosome block
scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) + # set y axis to 2 + log10 min p value that isn't zero
scale_color_manual(values = color_list)+
scale_size_continuous(range = c(0.5,3)) +
labs(x = NULL,
     y = "-log10(P)",color="Gene Type") +
theme_minimal() +
theme(
  #legend.position = "none",
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
)  + geom_label_repel(data=subset(mets_hits, is_annotate=="yes"), aes(label=Gene_label), size=4, force=10,show.legend = FALSE)+ # subset(mets_hits, is_annotate=="yes" | (Gene_Type_full!="Not Manually Annotated" & P_BOLT_LMM<=5e-8))
guides(color = guide_legend(override.aes = list(size=6)))

# Save plot
ggsave(paste0("manhattan_genetype_",group,".png"),width=12,height=8)}

met_groups=c("Ketone","BCAA","Amino_Acid","Glycolysis") # "Ketone","Amino_Acid","Glycolysis"

for (met_group in met_groups){
mets <- filter(cat,Color_Group==met_group)$Name_in_TSV
manh_plot(mets)
print(met_group)}
