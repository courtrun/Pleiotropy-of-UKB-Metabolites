# Courtney Smith - Nightingale Analysis - BOLT-LMM - Metabolite GWAS hit annotation
# Goal of script: Annotate closest gene + distance, distance to annotated gene minimum and TSS etc

#ml biology
#ml bedtools
#bedtools closest -a rigidmetabolites_sig_pruned_snplist.sorted.bed -b Homo_sapiens.GRCh37.87.chr.msiggenes.bed.gtf.gz -d > rigidmetabolites_sig_pruned_snplist.sorted.nearestgene.bed

library(dplyr)
g <- data.table::fread("Homo_sapiens.GRCh37.87.chr.genes.bed.gtf.gz")
an <- data.table::fread("rigidmetabolites_sig_pruned_snplist_manual_annotation.tsv")
ng <- data.table::fread("rigidmetabolites_sig_pruned_snplist.sorted.nearestgene.bed")

g <- g %>% mutate(tss=as.numeric(ifelse(strand=="+",chromStart,chromEnd)))
s <- tidyr::separate_rows(an,"Gene",sep=",") # into multiple rows
s <- left_join(s,g,by=c("Gene"="name"))
s <- s %>% mutate(Dist_Gene_TSS=abs(tss-BP),Within_Gene=ifelse(BP>chromEnd,"N",ifelse(BP<chromStart,"N","Y")))
s <- s %>% mutate(Dist_start=abs(BP-chromStart),Dist_end=abs(BP-chromEnd),Dist_min=pmin(Dist_start,Dist_end))
s <- s %>% mutate(Dist_Gene_Min=ifelse(Within_Gene=="Y",0,Dist_min))
s <- s %>% select(CHROM,BP,SNP,Variant_Type,Gene,Gene_Type,Dist_Gene_TSS,Dist_Gene_Min)
ng <- distinct(ng %>% select(V4,V8,V12) %>% group_by(V4) %>% arrange(V8)) %>% summarize(Nearest_Gene=paste(V8,collapse=","),Dist_Nearest_Gene=paste(V12,collapse=","))
dot <- s %>% filter(is.na(Dist_Gene_TSS))
s <- as.data.frame(s %>% group_by(SNP) %>% filter(Dist_Gene_Min == min(Dist_Gene_Min))%>% filter(1:n() == 1))
s <- bind_rows(s,dot)
s <- distinct(as.data.frame(distinct(s %>% group_by(SNP)) %>% arrange(Gene) %>% summarize(CHROM,BP,Gene=paste(Gene,collapse=","),Gene_Type,Variant_Type,Dist_Gene_TSS=paste(Dist_Gene_TSS,collapse=","),Dist_Gene_Min=paste(Dist_Gene_Min,collapse=","))))
s <- left_join(s,ng,by=c("SNP"="V4"))
s <- distinct(as.data.frame(distinct(s %>% group_by(SNP)) %>% summarize(CHROM,BP,Variant_Type,Gene,Gene_Type,Dist_Gene_TSS=min(Dist_Gene_TSS),Dist_Gene_Min,Nearest_Gene,Dist_Nearest_Gene)))

sfinal <- s %>% rename(Variant_Classification=Variant_Type,Assigned_Gene=Gene,Assigned_Gene_Type=Gene_Type,Dist_Assigned_Gene_TSS=Dist_Gene_TSS,Dist_Assigned_Gene_Min=Dist_Gene_Min)

OUTPUT_FILE="rigidmetabolites_sig_pruned_snplist_manual_annotation_closestgene_disttss.tsv"
write.table(sfinal, OUTPUT_FILE, quote=F, sep="\t", row.names=F, col.names=T)

(nrow(sfinal %>% filter(Nearest_Gene==Assigned_Gene))+6)/nrow(sfinal) # 69 %, added in the ZNF one that is missing from gene file and the genes where assigned was = to one of the multiple nearest genes assigned to that variant
# Dist_Nearest_Gene is from bedtools closet which uses least genomic distance from start or end of gene boundaries
