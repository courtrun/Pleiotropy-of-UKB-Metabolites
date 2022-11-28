
module load plink2

cat PCCB_VOI.txt | while read snp gene id; do
    mkdir -p conditional/${gene}clinical/
    grep -F $gene msig_within100kbOfGenes.all.bed
    grep -F $gene msig_within100kbOfGenes.all.bed | while read chrom start end name; do
        echo chromosome $chrom position $start to $end 
        test -f conditional/${gene}clinical/none.log || plink2 --glm cols=chrom,pos,ref,alt,a1freq,firth,test,nobs,orbeta,se,ci,tz,p hide-covar omit-ref --pfile ukbb_genotypes/ukb_imp_${chrom}_v3.mac1 --covar covar.phe --keep <(cut -f 1-2 Residual_extra.tsv) --out conditional/${gene}clinical/none --pheno Residual_extra.tsv --extract <(zstdcat ukbb_genotypes/ukb_imp_${chrom}_v3.mac1.pvar | awk "\$2 > $start - 1e6 && \$2 < $end + 1e6" | cut -f 3) --maf 0.01
        test -f conditional/${gene}clinical/${snp}.log || plink2 --glm cols=chrom,pos,ref,alt,a1freq,firth,test,nobs,orbeta,se,ci,tz,p hide-covar omit-ref --pfile ukbb_genotypes/ukb_imp_${chrom}_v3.mac1 --covar covar.phe --keep <(cut -f 1-2 Residual_extra.tsv) --out conditional/${gene}clinical/$snp --pheno Residual_extra.tsv --extract <(zstdcat ukbb_genotypes/ukb_imp_${chrom}_v3.mac1.pvar | awk "\$2 > $start - 1e6 && \$2 < $end + 1e6" | cut -f 3) --condition $id --maf 0.01
    done
done

