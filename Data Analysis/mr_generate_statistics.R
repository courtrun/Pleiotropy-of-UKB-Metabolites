
# Our pipeline is largely based on 
# https://academic.oup.com/ije/article/47/4/1264/5046668, Box 3
# Aka the Rücker model-selection framework
#
# This script is derived from the MR pipeline developed as part of
# Sinnott-Armstrong*, Tanigawa* et al. Nature Genetics, and was first
# written by Nasa Sinnott-Armstrong with help from David Amar. This
# version was updated to analyse the Nightingale metabolite data
# together with Courtney Smith.

library(TwoSampleMR)
library(data.table)
library(dplyr)

read_in_exposure <- function (d) {
    read_exposure_data(d, snp_col="SNP",
                          beta_col="BETA",
                          se_col="SE",
                          eaf_col="A1FREQ",
                          effect_allele_col="ALLELE1",
                          other_allele_col="ALLELE0",
                          pval_col="P_BOLT_LMM",
                          id_col="SNP", sep="\t", clump=F) %>% mutate(exposure=d, SNP=id.exposure, id.exposure=d)
}


hit_list = commandArgs(TRUE)[4:length(commandArgs(TRUE))]

exposure_dat <- rbindlist(lapply(hit_list, read_in_exposure))

outcome <- commandArgs(TRUE)[3]
outcome_raw <- data.table::fread(outcome, sep="\t", header=T, stringsAsFactors=F) %>% na.omit

disease.key <- read.table("mr_rg_traits.txt", header=F, col.names=c("name", "path", "title", "pmid"), sep="\t")

disease.key$sign = 1
disease.key$sign[grepl("ukb", ignore.case=T, disease.key$path)] <- -1
disease.key$sign[grepl("irnt", ignore.case=T, disease.key$path)] <- 1
disease.key$sign[grepl("X_RESPIRATORY", ignore.case=T, disease.key$path)] <- 1


print(head(exposure_dat$SNP))
print(head(outcome_raw$SNP))

#all_traits <- data.table::fread("../../paper/main/lcv_mr_rg_traits.txt", header=F, stringsAsFactors=F)
outcome_raw <- outcome_raw %>% filter(SNP %in% exposure_dat$SNP)

if (mean(outcome_raw$SIGNED_SUMSTAT, na.rm=T) > 0.5 && mean(outcome_raw$SIGNED_SUMSTAT < 0, na.rm=T) < 0.05) {
# convert ORs to betas on traits where the mean is greater than 0.5 and less than 5% are negative
    outcome_raw$SIGNED_SUMSTAT <- log(outcome_raw$SIGNED_SUMSTAT)
}

print(head(outcome_raw))

outcome_raw$P = 2*pnorm(-abs(outcome_raw$Z))

outcome_raw$SE=outcome_raw$SIGNED_SUMSTAT/outcome_raw$Z


outcome_dat <- format_data(as.data.frame(outcome_raw), type = "outcome", snps=unique(exposure_dat$SNP), snp_col="SNP", beta_col="SIGNED_SUMSTAT", se_col="SE", effect_allele_col="A1", other_allele_col="A2", samplesize_col="N", pval_col="P", min_pval=1e-200)
print(head(outcome_dat))
print(head(exposure_dat))

dat <- harmonise_data(exposure_dat, outcome_dat)

twosample <- c("mr_egger_regression", "mr_egger_regression_bootstrap", "mr_ivw_mre", "mr_ivw", "mr_ivw_fe")#, "mr_raps")

print("Saving data")
print(commandArgs(TRUE)[1])

save(dat, file=commandArgs(TRUE)[1])

res <- mr(dat, method_list=twosample, parameters=list(nboot=1000000))

print("Replacing with data + results")
save(res, dat, file=commandArgs(TRUE)[1])

pleiotropy <- mr_pleiotropy_test(dat)

print("Replacing with data + results + pleiotropy")
save(res, dat, pleiotropy, file=commandArgs(TRUE)[1])

steiger <- directionality_test(dat)

print("Replacing with data + results + pleiotropy + direction")
save(res, dat, pleiotropy, steiger, file=commandArgs(TRUE)[1])

heterogeneity <- mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))

print("Replacing with data + results + pleiotropy + direction + heterogeneity")
save(res, dat, pleiotropy, steiger, heterogeneity, file=commandArgs(TRUE)[1])
