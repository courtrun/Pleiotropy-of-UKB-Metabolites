# Courtney Smith - Nightingale Analysis - BOLT-LMM - Beta vs beta
# Goal of script: Plot the betas of one metabolite vs the betas of another for all nightingale sig pruned variants

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

options(ggrepel.max.overlaps = Inf)

rigidmet <- data.table::fread("mGWAS_rigidmetabolites_sig_pruned.tsv")
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation.tsv")

gencor <- data.table::fread("combined_gen_cor_sumstats.txt",fill=TRUE)
gencor <- na.omit(gencor) # get rid of pairs where one or more files did not exist
# Reformat trait columns to just be metabolite names
gencor$p1 <- gsub(".*Pass.","",gencor$p1)
gencor$p1 <- gsub(".txt.*","",gencor$p1)
gencor$p2 <- gsub(".*Pass.","",gencor$p2)
gencor$p2 <- gsub(".txt.*","",gencor$p2)

# Drop pairs of metabolites with themselves
gencor <- gencor %>% filter(p1!=p2) %>% select(p1,p2,rg,se,p,z)

######### Rigid mets
mets=c("Ala","Ile")
matrix_name="rigidmet"
snp_to_label1 = "rs370014171"
Gene_name1 = ": PDPR"
snp_to_label2 = "rs77010315"
Gene_name2 = ": SLC36A2"

met1=mets[1]
met2=mets[2]
beta1 <- paste0("BETA_",met1)
beta2 <- paste0("BETA_",met2)
pval1 <- paste0("P_BOLT_LMM_",met1)
pval2 <- paste0("P_BOLT_LMM_",met2)
se1 <- paste0("SE_",met1)
se2 <- paste0("SE_",met2)
p_partial_thres = 1e-4
p_sig_thres = 5e-8
matrix=rigidmet

keep <- matrix %>% select(matches(met1),matches(met2),SNP) %>%
				rename(BETA_met1=!!sym(beta1),BETA_met2=!!sym(beta2),P_BOLT_LMM_met1=!!sym(pval1),P_BOLT_LMM_met2=!!sym(pval2),SE_met1=!!sym(se1),SE_met2=!!sym(se2)) %>%
				mutate(p_partial_thres = ifelse(P_BOLT_LMM_met1<p_partial_thres & P_BOLT_LMM_met2<p_partial_thres,"Yes","No"),beta_opp=ifelse(sign(BETA_met1)==sign(BETA_met2),"No","Yes"))%>%
				mutate(p_sig_thres = ifelse((P_BOLT_LMM_met1<p_sig_thres | P_BOLT_LMM_met2<p_sig_thres)&(p_partial_thres=="Yes"),"Yes","No"))

# Plotting axes limits
bound <- max(abs(c(keep$BETA_met1,keep$BETA_met2)))+0.00001
x_down <- -0.037 # min(keep$met1_down)-0.0005
x_up <- 0.0262 # max(keep$met1_up)+0.0005
y_down <- -0.0282 # min(keep$met2_down)-0.0005
y_up <- 0.0365 # max(keep$met2_up)+0.0005

filename <- paste0(met1,"_vs_",met2,"_multiSNPlabel")
gc <- filter(gencor,p1==met1&p2==met2)$rg
keep <- keep %>% mutate(met1_down=BETA_met1-SE_met1,met1_up=BETA_met1+SE_met1,met2_down=BETA_met2-SE_met2,met2_up=BETA_met2+SE_met2)
keep <- keep %>% mutate(SNP_label=ifelse(SNP %in% c(snp_to_label1,snp_to_label2),SNP,""))
keep <- keep %>% mutate(Gene_label=ifelse(SNP==snp_to_label1,Gene_name1,ifelse(SNP==snp_to_label2,Gene_name2,"")),beta_opp_sign=ifelse(beta_opp=="Yes",-1,1))

ggplot(keep,aes(x=BETA_met2,y=BETA_met1))+
	geom_point(aes(alpha=p_partial_thres,color=ifelse(sign(beta_opp_sign)!=sign(gc) & p_partial_thres=="Yes","Yes","No")))+
	labs(y="Beta of Alanine",x="Beta of Isoleucine",color="Discordant Variant?",alpha="P < 1e-4 in both?")+
  scale_color_manual(values = c("black", "#5eba5e"))+ # #6bae6b #85e085
	ylim(x_down,x_up)+xlim(y_down,y_up) + # xlim(-bound,bound)+ylim(-bound,bound)+
  geom_errorbar(aes(ymin=met1_down, ymax=met1_up,alpha=p_partial_thres,color=ifelse(sign(beta_opp_sign)!=sign(gc) & p_partial_thres=="Yes","Yes","No"))) + # ,position=position_dodge(.9)
	geom_errorbarh(aes(xmin=met2_down, xmax=met2_up,alpha=p_partial_thres,color=ifelse(sign(beta_opp_sign)!=sign(gc) & p_partial_thres=="Yes","Yes","No"))) + # ,position=position_dodge(.9)
  geom_hline(yintercept=0, linetype='dotted')+
	geom_vline(xintercept=0, linetype='dotted')+
  annotate("text", label = paste("LDSC Genetic correlation:\n",gc),x = -bound/2, y = bound-0.005, size = 4, color = "black")+
	geom_smooth(method='lm',
              color="black")+coord_fixed()+
	geom_label_repel(aes(label=paste0(SNP_label,Gene_label),color=ifelse(sign(beta_opp_sign)!=sign(gc) & p_partial_thres=="Yes","Yes","No")), size=4, min.segment.length=unit(0,'lines'), show.legend = FALSE,force=10,segment.size=0.25)+ # ,nudge_x=max/4,nudge_y=max/2, data=filter(loads_path,Lipoprotein!="N")
	theme_classic()+ theme(text = element_text(size = 20))
print(ggsave(paste0(matrix_name,filename,".pdf")))
