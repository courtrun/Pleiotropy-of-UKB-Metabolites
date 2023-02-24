# This script takes the combined pathway results and calculates the total genetic covariance within pathways

# pathwayfile si the combined pathway lists that give the genomic coordinates of every pathway
pathway_genes_dir="/"
pathwayfile=${pathway_genes_dir}/combined_pathway_wlipids_ureamet_tfandtransporters.txt

# load bedtools if needed on your system
module load bedtools

# output files
rm pathway_group_rhog_tf_transporter_minmaf05.txt
rm pathway_group_h2_tf_transporter_minmaf05.txt

# pathway group specific genetic covariance
for pg in `grep -v Pathway_Group $pathwayfile | cut -f 6 | sort | uniq`; do
    for met in local_rhog/step3_minmaf05/*.h2_snps.bed; do
        intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep $pg )  -wa | sortBed -i | uniq | awk -v met=`basename $met _step3.h2_snps.bed` -v pg=$pg '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
    done
done | tr ' ' '\t' >> pathway_group_rhog_tf_transporter_minmaf05.txt

# genome wide covariance
for met in local_rhog/step3_minmaf05/*.h2_snps.bed; do
    awk -v met=`basename $met _step3.h2_snps.bed` -v pg=GENOME '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}' $met
done | tr ' ' '\t' >> pathway_group_rhog_tf_transporter_minmaf05.txt

# Full combined pathway (minus transporters / transcription factors / lipids) covariance
for met in local_rhog/step3_minmaf05/*.h2_snps.bed; do
    intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep -v TRANSMEMBRANE_TRANSPORT | grep -v TRANSCRIPTION | grep -v LIPID )  -wa | sortBed -i | uniq | awk -v met=`basename $met _step3.h2_snps.bed` -v pg=ALL '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
done | tr ' ' '\t' >> pathway_group_rhog_tf_transporter_minmaf05.txt

# detailed pathway specific covariance
for pg in `grep -v Pathway_Group $pathwayfile | cut -f 4 | sort | uniq`; do
    for met in local_rhog/step3_minmaf05/*.h2_snps.bed; do
        intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep $pg )  -wa | sortBed -i | uniq | awk -v met=`basename $met _step3.h2_snps.bed` -v pg=$pg '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
    done
done | tr ' ' '\t' >> pathway_group_rhog_tf_transporter_minmaf05.txt

# outside-of-full-combined-pathway covariance
for met in local_rhog/step3_minmaf05/*.h2_snps.bed; do
    intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep -v TRANSMEMBRANE_TRANSPORT | grep -v TRANSCRIPTION | grep -v CUSTOM_TF | grep -v CUSTOM_TRANSPORTER | grep -v LIPID) -v -wa | sortBed -i | uniq | awk -v met=`basename $met _step3.h2_snps.bed` -v pg=NONE '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
done | tr ' ' '\t' >> pathway_group_rhog_tf_transporter_minmaf05.txt


# pathway group heritability
for pg in `grep -v Pathway_Group $pathwayfile | cut -f 6 | sort | uniq`; do
    for met in *.h2_snps.bed; do
        intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep $pg ) -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v met=`basename $met _step2.h2_snps.bed` -v pg=$pg '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
    done
done | tr ' ' '\t' >> pathway_group_h2_tf_transporter_minmaf05.txt

# genome wide heritability
for met in *.h2_snps.bed; do
    awk -v met=`basename $met _step2.h2_snps.bed` -v pg=GENOME '(NF == 6) {n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}' $met
done | tr ' ' '\t' >> pathway_group_h2_tf_transporter_minmaf05.txt

# full combined pathway (except transporter / TF / lipid) heritability
for met in *.h2_snps.bed; do
    intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep -v TRANSMEMBRANE_TRANSPORT | grep -v TRANSCRIPTION | grep -v LIPID )  -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v met=`basename $met _step2.h2_snps.bed` -v pg=ALL '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
done | tr ' ' '\t' >> pathway_group_h2_tf_transporter_minmaf05.txt

# detailed pathway specific heritability
for pg in `grep -v Pathway_Group $pathwayfile | cut -f 4 | sort | uniq`; do
    for met in *.h2_snps.bed; do
        intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep $pg )  -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v met=`basename $met _step2.h2_snps.bed` -v pg=$pg '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
    done
done | tr ' ' '\t' >> pathway_group_h2_tf_transporter_minmaf05.txt

# outside-of-full-combined-pathway heritability
for met in *.h2_snps.bed; do
    intersectBed -a $met -b <(grep -v Pathway_Group $pathwayfile | grep -v TRANSMEMBRANE_TRANSPORT | grep -v TRANSCRIPTION | grep -v CUSTOM_TF | grep -v CUSTOM_TRANSPORTER | grep -v LIPID) -v -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v met=`basename $met _step2.h2_snps.bed` -v pg=NONE '{n+=$5; h2+=$4; var += $6*$6; i+=1} END {print met, pg, n, h2, sqrt(var), i}'
done | tr ' ' '\t' >> pathway_group_h2_tf_transporter_minmaf05.txt
