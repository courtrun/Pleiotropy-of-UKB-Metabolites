# Courtney Smith - Nightingale Analysis - BOLT-LMM - QQplots
# Goal of script: Making qqplots from pvalues of metaboltie GWAS hits

library("ggplot2")
library("ggdendro")
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggbiplot)
library(forcats)
library(ggrepel)

options(ggrepel.max.overlaps = Inf)

an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv")

# Combine Assigned_General Cell function into one caetegory
an$Assigned_Gene_Type <- ifelse(an$Assigned_Gene_Type=="General_Cell_Function_Lipid","General_Cell_Function",an$Assigned_Gene_Type)
full <- data.table::fread("mGWAS_rigidmetabolites_sig_pruned.tsv")

betas_full <- full %>% select(matches("BETA"))
colnames(betas_full) <- gsub("BETA_","",colnames(betas_full))
se_full <- full %>% select(matches("SE_"))
z <- betas_full/se_full
z[rowMedians(as.matrix(z))<0,] <- z[rowMedians(as.matrix(z))<0,]*-1 # flip sign so SNPs have median zscore across metabolites thats positive
pvalues <- full %>% select(SNP,matches("P_BOLT"))
z_full <- as.data.frame(z)
rownames(z_full) <- full$SNP
z_full$SNP <- rownames(z_full)
z.long <- reshape2::melt(z_full,id=c("SNP"))

z.long_an <- left_join(z.long,an,by=c("SNP"))
z_full_an <- left_join(z_full,an,by=c("SNP"))
z_full_an <- left_join(z_full_an,z_full_an %>% group_by(Assigned_Gene_Type) %>% dplyr::count() %>% dplyr::rename(snps_in_Assigned_Genetype=n),by=c("Assigned_Gene_Type"))
p_an <- left_join(pvalues,an,by=c("SNP")) %>% select(matches("P_BOLT"),Assigned_Gene_Type,SNP)
p_long <- reshape2::melt(p_an,id=c("Assigned_Gene_Type","SNP"))
p_long$value <- ifelse(p_long$value==0,1e-300,p_long$value)

# Drop the most significant trait for each snp since it was conditioned on having at least one trait highly significant
drop <- p_long %>% group_by(SNP) %>% arrange(value) %>% slice_head(n = 1)
p_long_d  <- anti_join(p_long, drop, by = c("SNP","variable"))

pvals <- p_long %>% arrange(value) %>% mutate(quantile=1:nrow(p_long)/nrow(p_long))
pvals <- pvals %>% mutate(exp_pval = quantile,neglog10_gwas = -log10(value),neglog10_expected = -log10(exp_pval))

dark_colors <- RColorBrewer::brewer.pal(n=8,name="Dark2")

### Calculate lambdas
lambdas <- list()
Assigned_Gene_types = unique(pvals$Assigned_Gene_Type)
for (Assigned_Gene_type in Assigned_Gene_types) {
temp <- filter(pvals,Assigned_Gene_Type==Assigned_Gene_type)
chisq = qchisq(temp$value,1,lower.tail=FALSE);
lambdas[[Assigned_Gene_type]] <- median(chisq) / qchisq(0.5,1)
}

lambda_list <- bind_rows(lambdas)
t_lambda <- as.data.frame(t(lambda_list)) %>% dplyr::rename(lambda=V1) %>% mutate(lambda=round(lambda,2))
t_lambda$Assigned_Gene_Type <- colnames(lambda_list)

pvals <- left_join(pvals,t_lambda,by=c("Assigned_Gene_Type")) %>% mutate(color_label=paste0(Assigned_Gene_Type,"\n(Lambda: ",lambda,")\n"))
pvals <- left_join(pvals,an %>% select(SNP,Assigned_Gene),by=c("SNP"))
pvals$Assigned_Gene_to_label <- ifelse(pvals$quantile<0.005,pvals$Assigned_Gene,"")

ggplot(pvals, aes(x = neglog10_expected, y = neglog10_gwas, color = factor(color_label))) + geom_point()+
geom_abline(intercept = 0, slope = 1)+ labs(color="Assigned_Gene_Type")+
scale_color_manual(values = c(dark_colors[2],dark_colors[8],dark_colors[4],dark_colors[3],dark_colors[1],dark_colors[6]))+
geom_label_repel(aes(label=Assigned_Gene_to_label), size=4, force=10,show.legend = FALSE)+
 theme_classic()+labs(x="Expected -log10(P)",y="Observed -log10(P)",color="Gene Type")+facet_wrap(~Assigned_Gene_Type)
ggsave("qqplots_facet.png")
