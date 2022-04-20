# Courtney Smith - Nightingale Analysis - BOLT-LMM - Phenotype Correlation Matrix
# Goal of script: Plot a correlation matrix of the residualized metabolite levels

library(data.table)
library(dplyr)
library(corrplot)

r <- data.table::fread("Residual_NMR_biomarkers.tsv")[,-c(1:2)]

corr <- cor(r, method = "pearson", use = "complete.obs")

corr <- as.data.frame(as.table(corr)) #turn into a 3-column table
corr <- na.omit(corr) #remove the NA values from above
keep_rigid <- c("Gly","Ala","Leu","Val","Ile","Pyruvate","Tyr","His","Phe","Gln","Acetone","bOHbutyrate","Acetoacetate","Glucose","Lactate","Citrate")
corr <- filter(corr,Var1%in%keep_rigid & Var2%in%keep_rigid)
corr <- filter(corr,Var1!=Var2)
mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq") #turn corr back into matrix in order to plot with corrplot

# Use correlation between variables as distance for ordering
dd <- dist(mtx_corr)
hc <- hclust(dd)
mtx_corr<-mtx_corr[hc$order, hc$order]
mtx_corr[upper.tri(mtx_corr,diag=TRUE)]<- NA # half to NAs

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
mets <- data.frame(mets=rownames(mtx_corr))
m <- left_join(mets,m,by=c("mets"="Name_in_TSV"))
cols <- c("#ff0000","#ff0000","#ffaa00","#ffaa00","#ffaa00","#00aaff","#00aaff","#00aaff","#ffaa00","#ffaa00","#aa00ff","#aa00ff","#aa00ff","#ffaa00","#ff0000","#ff0000")

colnames(mtx_corr) <- m$Biomarker_name
rownames(mtx_corr) <- m$Biomarker_name

#plot correlations visually
pdf("phenotype_correlation_matrix.rigidmets.pdf",height=18,width=15)
corrplot(mtx_corr,type="lower",tl.col=cols,tl.cex=3,col.lim = c(-1,1),cl.cex = 3,col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
	    "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE",
	    "#4393C3", "#2166AC", "#053061")))(200),order="hclust",is.corr=FALSE, method = c("square"),na.label=" ")
dev.off()
