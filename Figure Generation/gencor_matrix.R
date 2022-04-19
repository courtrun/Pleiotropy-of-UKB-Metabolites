# Courtney Smith - Nightingale Analysis - BOLT-LMM - LDSC Genetic Correlation
# Goal of script: Plot LDSC genetic correlation matrix for rigid mets x rigid mets; rigid mets x (all - rigid mets)

# ml R/4.0
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(corrplot)

combined_genetic_correlations="combined_gen_cor_sumstats.txt"
gencor <- data.table::fread(combined_genetic_correlations,fill=TRUE)

# get rid of pairs where one or more files did not exist
gencor <- na.omit(gencor)

# Reformat trait columns to just be metabolite names
gencor$p1 <- gsub(".*Pass.","",gencor$p1)
gencor$p1 <- gsub(".txt.*","",gencor$p1)
gencor$p2 <- gsub(".*Pass.","",gencor$p2)
gencor$p2 <- gsub(".txt.*","",gencor$p2)

# Switch this into wide matrix form
x2 = dcast(gencor , p1 ~ p2, value.var = "rg")
mat2 = as.matrix(x2[, 2:ncol(x2)])
rownames(mat2) = x2$p1
mat2[mat2 > 1] = 1
mat2[mat2 < -1] = -1
corr <- mat2

# Use correlation between variables as distance for ordering
dd <- dist(mat2)
hc <- hclust(dd)
mat2<-mat2[hc$order, hc$order]
mat2[upper.tri(mat2,diag=TRUE)]<- NA # half to NAs

cm <- as.data.frame(mat2) %>% mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>%
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE)
cm <- cm %>% rename(Metabolite1=Var1,Metabolite2=Var2)

corr <- as.data.frame(as.table(corr))
corr <- na.omit(corr) #remove the NA values
#sort by highest correlation
corr <- corr[order(-abs(corr$Freq)),]
mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
mtx_corr[is.na(mtx_corr)] <- 0

## Repeat above but only for correlations rigid Metabolites
keep_rigid <- c("Gly","Ala","Leu","Val","Ile","Pyruvate","Tyr","His","Phe","Gln","Acetone","bOHbutyrate","Acetoacetate","Glucose","Lactate","Citrate")
corr_rigid <- filter(corr,Var1%in%keep_rigid & Var2%in%keep_rigid)
corr_rigid <- filter(corr_rigid,!(Var1==Var2))
mtx_corr <- reshape2::acast(corr_rigid, Var1~Var2, value.var="Freq")
mtx_corr[is.na(mtx_corr)] <- 0
cols=c("blue","red","green")

# Add coloring to metabolite labels
cat <- data.table::fread("biomarkers_pathwayorder.csv") %>% select(Name_in_TSV,Color_Group,Biomarker_name)
cat <- as.data.frame(cat)
# Add BCAA category
BCAA_list <- c("Leu","Ile","Val")
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[1]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[2]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[3]] <- "BCAA"
m <- data.frame(Biomarker_name=c("3-Hydroxybutyrate","Acetoacetate","Acetone","Citrate","Lactate","Pyruvate","Phenylalanine","Tyrosine","Glucose","Alanine","Leucine","Isoleucine","Valine","Glycine","Glutamine","Histidine"))
m <- m %>% left_join(cat,by=c("Biomarker_name")) %>% mutate(met_color_label=ifelse(Color_Group=="Amino_Acid","#ffaa00",ifelse(Color_Group=="BCAA","#00aaff",ifelse(Color_Group=="Glycolysis","#ff0000","#aa00ff"))))
cols <- m$met_color_label

colnames(mtx_corr) <- (m %>% arrange(Name_in_TSV))$Biomarker_name
rownames(mtx_corr) <- (m %>% arrange(Name_in_TSV))$Biomarker_name

#plot correlations visually
pdf("LDSC.gencorrplot.corrplot.rigidmets.pdf",height=18,width=15)
corrplot(mtx_corr,type="lower",tl.col=cols,tl.cex=3,col.lim = c(-1,1),cl.cex = 3,col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
	    "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE",
	    "#4393C3", "#2166AC", "#053061")))(200),order="hclust",is.corr=FALSE, method = c("square"),na.label=" ")
dev.off()
