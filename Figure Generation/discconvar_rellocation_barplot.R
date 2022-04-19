# Courtney Smith - Nightingale Analysis - BOLT-LMM - Discordant variant analysis
# Goal of script: Bar plot of relative location by variant type enrichment

library(dplyr)

an <- data.table::fread("discorvar_rellocation.tsv")

### Run analysis
an <- an %>% left_join(an %>% group_by(Variant_Classification,Relative_Location) %>% count(),by=c("Variant_Classification","Relative_Location"))
an <- an %>% left_join(an %>% group_by(Variant_Classification) %>% count(),by=c("Variant_Classification")) %>% rename(num_invartyperelloc=n.x,num_in_vartype=n.y)
an <- an %>% mutate(frac_per_vartype=num_invartyperelloc/num_in_vartype)
an <- an %>% mutate(se_vt = sqrt(num_invartyperelloc*(num_in_vartype-num_invartyperelloc)/(num_in_vartype^3)))#
an$Relative_Location <- as.factor(an$Relative_Location)
an$Variant_Classification <- as.factor(an$Variant_Classification)
an <- distinct(an %>% select(-c("SNP","Assigned_Gene","Assigned_Gene_Type","Metabolite.1","Metabolite.2")))

library("ggplot2")

ggplot(an,aes(y=Variant_Classification,x=frac_per_vartype,fill=Relative_Location))+ geom_bar(stat="identity", position=position_dodge2(reverse = TRUE),width=0.6)+
	labs(y="Variant Type",x="Fraction of Variants",fill="Relative Location")+
scale_fill_manual(values=c("#00A843","black"))+ # #85e085,
	geom_errorbar(aes(xmin=frac_per_vartype-se_vt*2, xmax=frac_per_vartype+se_vt*2,y=Variant_Classification),  width=.6,position=position_dodge2(width = .9,reverse = TRUE, padding = .6))+
	theme_classic()+theme(text=element_text(size=16)) # axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
ggsave("barplot_disconrellocation.pdf",height=4,width=9)

# for dis and con vars annotated with pathway relevant enzyme genes:
# num concordant with all met pairs "not between" = 12, num concordant with at least one met pair "between" = 2, num discordant with at least one met pair "not between" = 1, num discordant with all met pairs "between" = 5
fisher.test(matrix(c(12,2,1,5), nrow=2)) # 22.96 odds ration P = 0.0072
fisher.test(matrix(c(5,1,2,12), nrow=2)) # 22.96 odds ration P = 0.0072
fisher.test(matrix(c(2,12,5,1), nrow=2)) # # 0.04355 odds ratio, P = 0.007224
