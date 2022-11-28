Folder for scripts analyzing the processed GWAS summary statistics

gene_metgroup_assignments.R: Assigns each gene annotation of a metabolite GWAS hit to its biochemical group

pathway_regions_localrho.R: Makes pathway regions for local genetic correlation analysis

MR Scripts
==========
Our pipeline is largely based on the Rücker model-selection framework https://academic.oup.com/ije/article/47/4/1264/5046668, Box 3.

The files needed to repeat these analyses are shown below:

mr_generate_statistics.R: run MR and export all regression varieties needed for the Rucker model selection framework.

mr_export_betas.R: generate the tabular results for betas from MR regressions.

mr_export_heterogeneity.R: generate the tabular results for heterogeneity statistics from MR regressions.

mr_export_intercepts.R: generate the tabular results for intercepts from MR regressions.

mr_result_processing.R: apply the Rucker model selection framework to identify which MR regression statistics to report.

Conditional Analysis Scripts
============================
There are a number of scripts used for the conditional analyses. These scripts and their associated input files are located in the Conditional folder.
