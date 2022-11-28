
module load plink2

cat VOI_within500kb.txt | while read snp gene id; do
    mkdir -p conditional/$gene/
    grep -F $gene msig_within100kbOfGenes.all.bed
    grep -F $gene msig_within100kbOfGenes.all.bed | while read chrom start end name; do
        echo chromosome $chrom position $start to $end 
        test -f conditional/$gene/none.log || plink2 --glm cols=chrom,pos,ref,alt,a1freq,firth,test,nobs,orbeta,se,ci,tz,p hide-covar omit-ref --pfile ukbb_genotypes/ukb_imp_${chrom}_v3.mac1 --covar covar.phe --keep <(cut -f 1-2 ../technical_variation/Residual_NMR_biomarkers.tsv | grep -vFf relatives/european_phe_and_remove.phe) --out conditional/$gene/none --pheno ../greml/Residual_rigid_mets.tsv --extract <(zstdcat ukbb_genotypes/ukb_imp_${chrom}_v3.mac1.pvar | awk "\$2 > $start - 1e6 && \$2 < $end + 1e6" | cut -f 3)
        plink2 --glm cols=chrom,pos,ref,alt,a1freq,firth,test,nobs,orbeta,se,ci,tz,p hide-covar omit-ref --pfile ukbb_genotypes/ukb_imp_${chrom}_v3.mac1 --covar covar.phe --keep <(cut -f 1-2 ../technical_variation/Residual_NMR_biomarkers.tsv | grep -vFf relatives/european_phe_and_remove.phe) --out conditional/$gene/$snp --pheno ../greml/Residual_rigid_mets.tsv --extract <(zstdcat ukbb_genotypes/ukb_imp_${chrom}_v3.mac1.pvar | awk "\$2 > $start - 1e6 && \$2 < $end + 1e6" | cut -f 3) --condition $id
    done
done
