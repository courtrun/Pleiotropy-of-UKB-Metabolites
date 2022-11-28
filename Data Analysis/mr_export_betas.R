# Our pipeline is largely based on 
# https://academic.oup.com/ije/article/47/4/1264/5046668, Box 3
# Aka the Rücker model-selection framework
#
# This script is derived from the MR pipeline developed as part of
# Sinnott-Armstrong*, Tanigawa* et al. Nature Genetics, and was first
# written by Nasa Sinnott-Armstrong with help from David Amar. This
# version was updated to analyse the Nightingale metabolite data
# together with Courtney Smith.

library(dplyr)

mrfiles <- as.character(gsub(".RData", "", grep("RData", list.files("mrbase_local/"), value=T)))

disease.key <- read.table("mr_rg_traits.txt", header=F, col.names=c("name", "path", "title", "pmid"), sep="\t")

all.res <- data.table::rbindlist(lapply(mrfiles, function(x) {load(sprintf("mrbase_local/%s.RData", x)); if (!is.null(res) && !is.na(res) && nrow(res) > 0) {res$outcome = x; res}else{NULL}})) %>% 
    inner_join(disease.key, by=c("outcome" = "name"))

all.res$exposure <- gsub("[.]hits$", "", all.res$exposure)
all.res$exposure <- gsub("[.][.].*hits.", "", all.res$exposure)

# rescale betas for dichotomous traits (not applicable)
all.res$scaled.b <- ifelse(!is.na(all.res$ncase), log((all.res$ncase/all.res$ntotal + all.res$b)/(all.res$ncase/all.res$ntotal)), all.res$b)

all.res %>%
    write.table("all.results.restructured.txt", quote=F, sep="\t", row.names=F, col.names=T)
