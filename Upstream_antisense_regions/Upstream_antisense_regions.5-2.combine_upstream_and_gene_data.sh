#!/bin/bash

CASE=$1
CTRL=$2
FEATUREFOLDER=$3
BAMFOLDER=$4
NEIGHBORHOODFOLDER=$5

CTRLRPM=$(samtools idxstats $BAMFOLDER/$CTRL.best.bam | awk '{sum+=$3}END{print sum/1000000}')
CASERPM=$(samtools idxstats $BAMFOLDER/$CASE.best.bam | awk '{sum+=$3}END{print sum/1000000}')

join -t$'\t' -a 1 -e "#N/A" -o 1.1,2.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.3 <(awk -v rpm=$CTRLRPM 'NR==1{print $1"\t"$2"\t"$3"\t"$4"\tctrl_"$5"\tctrl_"$6"\tctrl_"$7"\tctrl_"$8"\tctrl_"$9"\tctrl_"$10"\tctrl_"$11"\tctrl_"$12}NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5/rpm"\t"$6"\t"$7/rpm"\t"$8/rpm"\t"$9/rpm"\t"$10/rpm"\t"$11/rpm"\t"$12/rpm}' $FEATUREFOLDER/$CTRL.feature) <(sort -k1,1 $CTRL.upstream_norm_reads.txt | sed 's|_upstream||' | cat <(echo -e "gene_ID\tgene_name\tctrl_upstream_antisense_reads") -) > $CTRL.feature_info.txt

join -t$'\t' -a 1 -e "#N/A" -o 1.1,2.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.3 <(awk -v rpm=$CASERPM 'NR==1{print $1"\t"$2"\t"$3"\t"$4"\tcase_"$5"\tcase_"$6"\tcase_"$7"\tcase_"$8"\tcase_"$9"\tcase_"$10"\tcase_"$11"\tcase_"$12}NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5/rpm"\t"$6"\t"$7/rpm"\t"$8/rpm"\t"$9/rpm"\t"$10/rpm"\t"$11/rpm"\t"$12/rpm}' $FEATUREFOLDER/$CASE.feature) <(sort -k1,1 $CASE.upstream_norm_reads.txt | sed 's|_upstream||' | cat <(echo -e "gene_ID\tgene_name\tcase_upstream_antisense_reads") -) > $CASE.feature_info.txt

paste $CTRL.feature_info.txt <(cut -f5- $CASE.feature_info.txt) > $CASE.tmp.txt

python /home/dreyfusslab/Desktop/projects/gene_neighborhoods/upstream_antisense/combine_ctrl_case_upstream_antisense.py $NEIGHBORHOODFOLDER/with_intronic_pcpa/intermediate/$CASE.combined_pcpaed_genes.IDs $NEIGHBORHOODFOLDER/intermediate/$CASE.no_exon1_upregulated.IDs $NEIGHBORHOODFOLDER/with_intronic_pcpa/intermediate/$CASE.down_not_PCPAed.IDs $CASE.tmp.txt > $CASE.upstream_antisense.txt

rm $CASE.tmp.txt
