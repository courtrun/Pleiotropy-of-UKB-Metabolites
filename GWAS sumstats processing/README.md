Folder for scripts processing GWAS sumstats

Snakefile_gwasprocessing: BOLT-LMM nightingale metabolite GWAS sumstats processing Snakefile; calls below R scripts

mfs.R: Filters by MAF and INFO; Called in Snakefile_gwasprocessing by rule filter_plink

filtersumstats_toprunedsnps.R: Filter sumstats to just pruned snplist; called in Snakefile_gwasprocessing script by rule pruned_snplist

make_matrix.R: Make matrix of pruned hits by rigid metabolites; Called in Snakefile_gwasprocessing by rule make_matrix
