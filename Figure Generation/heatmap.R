# Courtney Smith - Nightingale Analysis - BOLT-LMM - Rigid Metabolite Heatmap Plot
# Goal of script: Plot heatmap of zscores for rigid metabolites vs variants ascertained in them; with and without x axis

library("ggplot2")
library("ggdendro")
library("reshape2")
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggbiplot)
library(forcats)

full <- data.table::fread("mGWAS_rigidmetabolites_sig_pruned.tsv")

betas_full <- full %>% select(matches("BETA"))
colnames(betas_full) <- gsub("BETA_","",colnames(betas_full))
se_full <- full %>% select(matches("SE_"))
z <- betas_full/se_full
z <- t(z)
my_function <- function(x) {return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}
z <- t(apply(z,FUN=my_function,MARGIN=1)) # apply transform to each row (1=row, 2=column), row normalize
means <- rowMeans(as.matrix(z)) # calculate the mean of each row (aka each metabolite) of z
SDs <- rowSds(as.matrix(z)) # calculate the SD of each row (aka each metabolite) of z
z <- sweep(z, 1, means) # subtract each row by its respective mean
z <- sweep(z, 1, SDs, FUN = '/') # and divide each row by its respective SD
colnames(z) <- full$SNP
z <- t(z)
z[rowMedians(as.matrix(z))<0,] <- z[rowMedians(as.matrix(z))<0,]*-1 # flip sign so SNPs have median zscore across metabolites thats positive
pvalues <- full %>% select(SNP,matches("P_BOLT"))
z_full <- as.data.frame(z)
z_full$SNP <- rownames(z_full)
z.long <- melt(z, id = c("Metabolite"))

set.seed(17) # untransposed, order SNPs
z.dendro <- as.dendrogram(hclust(d = dist(x = z, method = "euclidean"), method = "ward.D2"))
dendro.plot <- ggdendrogram(data = z.dendro, rotate = TRUE) +
  theme(axis.text.y = element_text(size = 30))
z.order <- order.dendrogram(z.dendro)

set.seed(17) # transposed, order metabolites
library("dendextend")
z.dendro2 <- as.dendrogram(hclust(d = dist(x = t(z), method = "euclidean"), method = "ward.D2"))
z.order2 <- order.dendrogram(z.dendro2)
# Add coloring to metabolite labels
cat <- data.table::fread("biomarkers_pathwayorder.csv") %>% select(Name_in_TSV,Color_Group,Biomarker_name)
cat <- as.data.frame(cat)
# Add BCAA category
BCAA_list <- c("Leu","Ile","Val")
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[1]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[2]] <- "BCAA"
cat["Color_Group"][cat["Name_in_TSV"]==BCAA_list[3]] <- "BCAA"
ddata_x <- dendro_data(z.dendro2)

# Dendrogram
p2 <- ggplot(segment(ddata_x)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+ coord_flip()+
      #theme(axis.text.y = element_text(size = 30),axis.text.x=element_blank())+
      theme_dendro()+ theme(legend.position = "none")

# Text
labs <- label(ddata_x)
labs <- labs %>% left_join(cat,by=c("label"="Name_in_TSV")) %>% mutate(met_color_label=ifelse(Color_Group=="Amino_Acid","#ffaa00",ifelse(Color_Group=="BCAA","#00aaff",ifelse(Color_Group=="Glycolysis","#ff0000","#aa00ff"))))
dendro.text <- ggplot()+geom_label(data=label(ddata_x) %>% left_join(cat,by=c("label"="Name_in_TSV")), aes(label=Biomarker_name, x=x, y=y+2.6, color=as.factor(labs$Color_Group)),size = 6)+
            scale_color_manual(values = unique((labs %>% arrange(Color_Group))$met_color_label))+ coord_flip()+
            #theme(axis.text.y = element_text(size = 30),axis.text.x=element_blank())+
            theme_dendro()+
             theme(legend.position = "none")

pdf('dendrogram.pdf')
grid.newpage()
print(p2, vp = viewport(x = 0.6, y = 0.5, width = 0.6, height = 1.0))
print(dendro.text, vp = viewport(x = 0.2, y = 0.5, width = 0.35, height = 1))
dev.off()

sig_thres = 1e-6
sub_sig_thres = 1e-3

p <- reshape2::melt(pvalues, id.vars="SNP") %>% dplyr::rename(Metabolites=variable,P_BOLT_LMM=value)
p$Metabolites <- gsub("P_BOLT_LMM_","",p$Metabolites)
z.long <- left_join(z.long,p,by=c("Var1"="SNP","Var2"="Metabolites"))

# Order the levels according to their position in the cluster
z.long$Var1 <- factor(x = z.long$Var1,
                               levels = z_full$SNP[z.order],
                               ordered = TRUE)

z.long$Var2 <- factor(x = z.long$Var2,
                              levels = colnames(z_full)[z.order2],
                              ordered = TRUE)

## Adding coloring by gene type
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation.tsv")
# Combine General Cell function into one caetegory
an$Gene_Type <- ifelse(an$Gene_Type=="General_Cell_Function_Lipid","General_Cell_Function",an$Gene_Type)
an <- an %>% mutate(SNP_label=paste0(SNP," - ",Gene))
z.long_an <- left_join(z.long,an,by=c("Var1"="SNP"))
z_full_an <- left_join(z_full,an,by=c("SNP"))

z.long_an$SNP_label <- factor(x = z.long_an$SNP_label,
                               levels = z_full_an$SNP_label[z.order],
                               ordered = TRUE)

z.long_an$Var2 <- factor(x = z.long_an$Var2,
                              levels = colnames(z_full)[z.order2],
                              ordered = TRUE)
z_full_an$SNP_label <- factor(x = z_full_an$SNP_label,
                               levels = z_full_an$SNP_label[z.order],
                               ordered = TRUE)

z_full_an <- z_full_an %>% arrange(SNP_label)
dark_colors <- RColorBrewer::brewer.pal(n=8,name="Dark2")
color_list <- ifelse(z_full_an$Gene_Type=="Transporter",dark_colors[1],ifelse(z_full_an$Gene_Type=="General_Cell_Function",dark_colors[2],ifelse(z_full_an$Gene_Type=="Transcription",dark_colors[3],ifelse(z_full_an$Gene_Type=="Pathway_Relevant_Enzyme",dark_colors[4],dark_colors[6]))))
color_df <- data.frame(Gene_Type=unique(z_full_an$Gene_Type),Color=unique(color_list),x=1:length(unique(z_full_an$Gene_Type)),y=1:length(unique(z_full_an$Gene_Type)))

# Plot without x axis
heatmap.plot <- ggplot(data = z.long_an, aes(x = SNP_label, y = Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = ifelse(P_BOLT_LMM < sig_thres, "*", "")), vjust = 1) +
  geom_text(aes(label = ifelse(P_BOLT_LMM > sig_thres & P_BOLT_LMM < sub_sig_thres, ".", "")))+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
midpoint = 0, name="Normalized Z-scores")+theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size = 70),
        legend.position = "top",legend.key.size = unit(3, 'cm'),
         legend.text = element_text(size=40),
         legend.title = element_text(size=70))+
      labs(x="SNPs",y="Metabolites")
pdf('euclidean_wardd2_annotated_genetype_noxaxis.pdf',width=45,height=20)
print(heatmap.plot)
dev.off()
