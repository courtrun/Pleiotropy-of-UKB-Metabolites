# Our pipeline is largely based on 
# https://academic.oup.com/ije/article/47/4/1264/5046668, Box 3
# Aka the Rücker model-selection framework
#
# This script is derived from the MR pipeline developed as part of
# Sinnott-Armstrong*, Tanigawa* et al. Nature Genetics, and was first
# written by David Amar with help from Nasa Sinnott-Armstrong. This
# version was updated to analyse the Nightingale metabolite data
# together with Courtney Smith.

library(data.table)
library(dplyr)
# read the trait name info
trait_names = fread("mr_rg_traits.txt",
                    data.table = F,stringsAsFactors = F,
                    header= F, sep="\t")

rigid <- data.table::fread("metadata/rigidmets.txt", stringsAsFactors=F, header=F, sep="\t")
all <- data.table::fread("metadata/metabolite_list.txt", stringsAsFactors=F, header=F, sep="\t")

nondisease <- data.table::fread("mr_rg_traits_nondisease.txt", stringsAsFactors=F, header=F, sep="\t")

rnames = trait_names$V1
trait_names = trait_names$V3
names(trait_names) = rnames

# read the results 
simplify_name<-function(s,trait_names){
  if(is.element(s,rownames(trait_names))){
    return(trait_names[s])
  }
  arr = strsplit(s,split="/|,")
  return(gsub(" ", "_", gsub(".*metabolomeFirstPass.", "", arr[[1]][length(arr[[1]])])))
}

raw_results = fread("./all.results.restructured.txt",stringsAsFactors = F,data.table = F)
raw_results$exposure = sapply(raw_results$exposure,simplify_name,trait_names=trait_names)
raw_results$outcome = sapply(raw_results$outcome,simplify_name,trait_names=trait_names)

pleio_res = fread("./all.restructured.intercept.txt",stringsAsFactors = F,data.table = F)
rownames(pleio_res) = paste(pleio_res$exposure,pleio_res$outcome,sep=",")
pleio_res$exposure = sapply(pleio_res$exposure,simplify_name,trait_names=trait_names)
pleio_res$outcome = sapply(pleio_res$outcome,simplify_name,trait_names=trait_names)
het_res = fread("./all.restructured.heterogeneity.txt",stringsAsFactors = F,data.table = F)
het_res$exposure = sapply(het_res$exposure,simplify_name,trait_names=trait_names)
het_res$outcome = sapply(het_res$outcome,simplify_name,trait_names=trait_names)

print(head(raw_results))

het_res <- filter(het_res, exposure %in% rigid$V1, exposure != outcome, !(outcome %in% all$V1), !(outcome %in% grep("Ratio", outcome, value=T)), !(outcome %in% nondisease$V1))
pleio_res <- filter(pleio_res, exposure %in% rigid$V1, exposure != outcome, !(outcome %in% all$V1), !(outcome %in% grep("Ratio", outcome, value=T)), !(outcome %in% nondisease$V1))
raw_results <- filter(raw_results, exposure %in% rigid$V1, exposure != outcome, !(outcome %in% all$V1), !(outcome %in% grep("Ratio", outcome, value=T)), !(outcome %in% nondisease$V1))

print(dim(het_res))
print(dim(pleio_res))
print(dim(raw_results))
print(head(raw_results))

raw_results = exclude_pairs(raw_results)
pleio_res = exclude_pairs(pleio_res)
het_res = exclude_pairs(het_res)

# Define the thresolds
Q_p_thr = 0.01
FDR_level = 1
# p_het_thr = 0.01

# separate the main results by method
method2res = list()
for(method in unique(raw_results$method)){
  method2res[[method]] = raw_results[raw_results$method==method,]
}
names(method2res) = gsub("Inverse variance weighted","IVW",names(method2res))
names(method2res) = gsub("random effects","RE",names(method2res))
names(method2res) = gsub("fixed effects","FE",names(method2res))
method2res = lapply(method2res,
    function(x){rownames(x) = paste(x$exposure,x$outcome,sep=",");
                x$exposure = sapply(x$exposure,simplify_name,trait_names=trait_names);
                x$outcome = sapply(x$outcome,simplify_name,trait_names=trait_names);x})
# Add the FDR to the methods
for(mname in names(method2res)){
  q_values = p.adjust(method2res[[mname]]$pval,method="fdr")
  method2res[[mname]] = cbind(method2res[[mname]],q_values)
}

########################################################################
########################################################################
# Method comparison:
# plot the raw p-values
par(mfrow=c(2,3))
for(method in names(method2res)){
  mname = strsplit(method,split="\\(|\\)")[[1]]
  mname = paste(mname,collapse="\n")
  hist(method2res[[method]]$pval,
       main=mname,xlab="P-value",cex.main=1)
}
# Method similarity
shared_pairs = rownames(method2res[[1]])
for(method in names(method2res)){
  shared_pairs = intersect(shared_pairs,rownames(method2res[[method]]))
}
pval_mat = sapply(method2res,function(x,y)x[y,"pval"],y=shared_pairs)
pval_mat = cbind(pval_mat,pleio_res[shared_pairs,]$pval)
colnames(pval_mat)[ncol(pval_mat)] = "Egger_pleio_p"
scaled_b_mat = sapply(method2res,function(x,y)x[y,"scaled.b"],y=shared_pairs)
library(ggcorrplot)
pval_corrs = cor(pval_mat,method="spearman")
print(ggcorrplot(t(pval_corrs),lab=T,lab_size=2.5,hc.order = F) +
        ggtitle("P-value, spearman") +
        theme(plot.title = element_text(hjust = 0.5,size=20)))
# remove rows with NAs
scaled_b_mat = scaled_b_mat[!apply(is.na(scaled_b_mat),1,any),]
b_corrs = cor(scaled_b_mat,method="spearman")
print(ggcorrplot(t(b_corrs),lab=T,lab_size=2.5,hc.order = F) +
        ggtitle("Scaled beta, spearman") +
        theme(plot.title = element_text(hjust = 0.5,size=20)))
dev.off()
########################################################################

# Adjust for FDR - get the number of results per method
get_adjusted_results<-function(x,sig=0.01){
  return(x[p.adjust(x$pval,method="fdr")<sig,])
}
method2adjres = lapply(method2res,get_adjusted_results)
par(mar=c(5,10,5,5))
barplot(sapply(method2adjres,nrow),las=2,horiz = T,xlab="Number of pairs (0.01 FDR)")
dev.off()
########################################################################
# Scheme for selecting the models and their results based on 
# the Rücker model-selection framework

# In our case (comparisons above) we used Egger with bootstrap to increas power.
# We then adapt the analysis of the heterogeneity as follows:
# For inignificant IVW Q scores, use the beta and p-value from IVW+FE.
# For significant Q scores use the beta valuse use IVW+ME.
# If the difference between the Egger Q and the IVW+FE Q is significant, we use the 
# beta and p-value from Egger.

# IVW
ivw_fe_results = method2res$`IVW (FE)`
ivw_fe_q_results = het_res[het_res$method == "Inverse variance weighted",]
ivw_me_results = method2res$IVW
rownames(ivw_fe_q_results) = paste(ivw_fe_q_results$exposure,ivw_fe_q_results$outcome,sep=",")
rownames(ivw_fe_results) = paste(ivw_fe_results$exposure,ivw_fe_results$outcome,sep=",")
rownames(ivw_me_results) = paste(ivw_me_results$exposure,ivw_me_results$outcome,sep=",")
ivw_shared_pairs = intersect(rownames(ivw_fe_q_results),rownames(ivw_fe_results))
ivw_fe_results = ivw_fe_results[ivw_shared_pairs,]
ivw_fe_q_results = ivw_fe_q_results[ivw_shared_pairs,]
ivw_me_results = ivw_me_results[ivw_shared_pairs,]
is_ivw_fe_q_significant = ivw_fe_q_results$Q_pval < Q_p_thr
names(is_ivw_fe_q_significant) = rownames(ivw_fe_q_results)
ivw_merged_results = rbind(
  ivw_fe_results[!is_ivw_fe_q_significant,],
  ivw_me_results[is_ivw_fe_q_significant,]
)

# Egger
# make sure that Egger and EggerB are ordered correctly: (sanity check)
all(method2res$`MR Egger`[,1:4] == method2res$`MR Egger (bootstrap)`[,1:4])
egger_results = method2res$`MR Egger`
egger_b_inds = egger_results$nsnp > 30
egger_results[egger_b_inds,] = method2res$`MR Egger (bootstrap)`[egger_b_inds,]
egger_q_results = het_res[het_res$method == "MR Egger",]
rownames(egger_results) = paste(egger_results$exposure,egger_results$outcome,sep=",")
rownames(egger_q_results) = paste(egger_q_results$exposure,egger_q_results$outcome,sep=",")
# Methods' rownames do not perfectly fit but all Egger pairs are in the IVW pairs
egger_q_diffs = ivw_fe_q_results[rownames(egger_q_results),"Q"] - egger_q_results$Q
egger_q_diffs_pval = pchisq(egger_q_diffs,1,lower.tail = F)
table(egger_q_diffs_pval > Q_p_thr)
egger_pairs = rownames(egger_q_results)[egger_q_diffs_pval < Q_p_thr]
egger_pairs = egger_pairs[egger_q_results[egger_pairs,"Q_pval"] > 1e-100]
ivw_egger_merged_results = ivw_merged_results
ivw_egger_merged_results[egger_pairs,] = egger_results[egger_pairs,]

# Get the current significant results
selected_merged_results = p.adjust(ivw_egger_merged_results$pval,method="fdr") < FDR_level
table(selected_merged_results)

selected_merged_results = ivw_egger_merged_results[selected_merged_results,]
rownames(selected_merged_results)
selected_merged_results = selected_merged_results[,-c(1:2)]

unique_pairs = unique(selected_merged_results[,2:1])
rownames(unique_pairs) = NULL
write.table(unique_pairs,sep="\t",quote=F,row.names = F)

# Add alternative beta scores
selected_merged_results[["abs(b)"]] = abs(selected_merged_results$scaled.b)
selected_merged_results[["log(b^2)"]] = log(selected_merged_results$scaled.b^2,base=10)
v = as.numeric(selected_merged_results$b>0)
v[v==0]=-1
selected_merged_results[["Effect_sign"]] = v

# Node attr for cytoscape
allnodes = union(selected_merged_results[,1],selected_merged_results[,2])
m = cbind(allnodes,is.element(allnodes,set=selected_merged_results$outcome))
colnames(m) = c("node","is_outcome")
write.table(m,file="disease_rigidExposure_node_info_fdr1_Q0.05.txt"
            ,sep="\t",row.names = F,col.names = T,quote = F)

write.table(selected_merged_results,
            file="disease_rigidExposure_results_fdr1_Q0.05.txt"
            ,sep="\t",row.names = F,col.names = T,quote = F)
