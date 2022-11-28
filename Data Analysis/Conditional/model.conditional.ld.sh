
module load plink2

cat VOI_within500kb.txt | while read snp gene id; do
    mkdir -p conditional/$gene/
    grep -F $gene msig_within100kbOfGenes.all.bed
    grep -F $gene msig_within100kbOfGenes.all.bed | while read chrom start end name; do
        echo chromosome $chrom position $start to $end 
	test -f ${gene}_ld.bed || plink2 --make-bed --pfile ukbb_genotypes/ukb_imp_${chrom}_v3.mac1 --keep <(cut -f 1-2 ../technical_variation/Residual_NMR_biomarkers.tsv | grep -vFf relatives/european_phe_and_remove.phe) --out ${gene}_ld --extract <(zstdcat ~/oak/ukbb_gt/ukb_imp_${chrom}_v3.mac1.pvar | awk -v start=$start -v end=$end '$2 > start - 1e6 && $2 < end + 1e6' | cut -f 3) --maf 0.01
         test -f ${gene}_ld_est.ld.gz || plink --bfile ${gene}_ld --out ${gene}_ld_est --r square gz --freq
    done
done
