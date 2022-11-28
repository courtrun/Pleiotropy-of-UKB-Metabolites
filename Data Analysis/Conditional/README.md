Folder for scripts running conditional analyses of loci identified through GWAS

`both_clinical_PCCBconditional.sh` - conditional GWAS of multiple variants at PCCB for clinical outcomes

`condition_PCCB_extended.R` - conditional GWAS of single variants at PCCB for the extended set of metabolites (fatty acids, LDL, etc)

`condition_PCCB.R` - conditional GWAS of single variants at PCCB for the extended set of metabolites (fatty acids, LDL, etc)

`condition_PDPR.R` - conditional GWAS of single variants at PDPR for the core set of metabolites

`condition_SLC36A2.R` - conditional GWAS of single variants at SLC36A2 for the core set of metabolites

`model.conditional.inputdata.sh` - script to generate the pgen files containing the subset of variants used in the conditional analyses

`model.conditional.ld.sh` - script to generate the LD matrices for conditional effect estimation

`model.conditional.PCCBclinical.sh` - script to run conditional analyses for PCCB clinical traits

`model.conditional.PCCB.sh` - script to run conditional analyses for PCCB metabolite traits

`model.conditional.sh` - script to run conditional analyses for other loci metabolite traits

`PCCB_VOI.txt` - list of PCCB variants of interest to test condition on

`VOI_within500kb.txt` - list of other loci variants of interest for conditional analyses

`PDPR_condition_plink_runs.sh` - script to generate & run conditioning for PDPR (treated specially due to its secondary association)

`run.all.coloc.PDPR_conditional.R` - script to run coloc on PDPR conditioned on the secondary signal

`run.all.coloc.PDPR.R` - script to run coloc on PDPR without extra conditioning

`run.all.coloc.SLC36A2.R` - script to run coloc on SLC36A2


