#!/bin/bash

### this script collects all of the junction information from TopHat and collates it into 3 groups:
# 1) junctions that splice using annotated 5' and 3' splice sites between consecutive exons (i.e. across known introns)
# 2) junctions that splice using annotated 5' and 3' splice sites but skip at least one intervening exon (i.e. multiexon skipping)
# 3) all other spliced reads in the gene (i.e. exon-intron, intron-exon, intra-intron or UTR, etc.)

## this takes the following inputs
## sample list
samples=$1 
## junction directory
juncDir=$2
## feature file directory
ffDir=$3
## sample comparison file with mapped reads in millions (if you don't want to run a comparison of case to ctrl input this as "0")
compFile=$4
## gff3 of all introns (ALL ISOFORMS, ex. mm10_all_introns-isoforms.gff3)
intronGFF3=$5
## gff3 of transcripts (RELEVANT ISOFORMS, ex. mm10_mrna_longest.gff3)
mrnaGFF3=$6
## longest isoform data file (name, ID, gene length, mRNA length, exon number, ex. mm10.longest_isoform.txt)
geneInfo=$7
## python script for getting multiexon skipping
pythonScript=get_multi_skipping_junctions.py


while read line
	do
	sed 's|,|\t|g' $juncDir/$line.junctions.bed | tail -n+2 | awk 'BEGIN{OFS=FS="\t"}{if(NF>4){print $1,$2+$13+1,$2+$16,$4,$5,$6}else{print $0}}' | awk '{if($5>2){print $0}}' > $line.junctions.bed.clean

	awk 'FNR==NR{a[$1,$4,$5]++;next}!a[$1,$2,$3]' $intronGFF3 $line.junctions.bed.clean > $line.junctions.aberrant

	python $pythonScript $intronGFF3 $line.junctions.aberrant $line.multi_junctions

	awk 'FNR==NR{a[$1,$4,$5]++;next}a[$1,$2,$3]' $intronGFF3 $line.junctions.bed.clean > $line.junctions.annotated

	awk 'FNR==NR{a[$4]++;next}!a[$4]' $line.multi_junctions $line.junctions.aberrant > $line.incorrect_splicing

	intersectBed -a $line.junctions.annotated -b $mrnaGFF3 -wo | cut -f1-6,15 | cut -f1 -d";" | sed 's|ID=||' > $line.junctions.annotated.clean

	intersectBed -a $line.incorrect_splicing -b $mrnaGFF3 -wo | cut -f1-6,15 | cut -f1 -d";" | sed 's|ID=||' > $line.incorrect_splicing.clean

	awk '{count[$7]++;sum[$7]+=$5}END{print "gene_ID\tannotated_junction_number\tannotated_junction_reads";for(i in count){print i"\t"count[i]"\t"sum[i]}}' $line.junctions.annotated.clean > $line.junctions.annotated.final
	awk '{count[$7]++;sum[$7]+=$5}END{print "gene_ID\tincorrect_junction_number\tincorrect_junction_reads";for(i in count){print i"\t"count[i]"\t"sum[i]}}' $line.incorrect_splicing.clean > $line.incorrect_splicing.final
	awk '{count[$7]++;sum[$7]+=$5}END{print "gene_ID\tmultiexon_junction_number\tmultiexon_junction_reads";for(i in count){print i"\t"count[i]"\t"sum[i]}}' $line.multi_junctions > $line.multi_junctions.final

	join -t$'\t' -1 2 -2 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 -a 1 -e "0" <(cat <(echo -e "gene_name\tgene_ID\tgene_length\ttranscript_length\texon_number") <(sort -k2,2 $geneInfo)) <(sort -k1,1 $line.junctions.annotated.final) | join -t$'\t' -1 2 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3 -a 1 -e "0" - <(sort -k1,1 $line.incorrect_splicing.final) | join -t$'\t' -1 2 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3 -a 1 -e "0" - <(sort -k1,1 $line.multi_junctions.final) > $line.per_gene_junction_breakdown.txt

	rm $line.multi_junctions $line.multi_junctions.final $line.incorrect_splicing.clean $line.incorrect_splicing.final $line.junctions.annotated.clean $line.junctions.annotated.final $line.junctions.aberrant $line.junctions.annotated $line.incorrect_splicing $line.junctions.bed.clean
	done < $samples
##
if [ $1 = 0 ]
    then
        echo "Skipping case vs. control analysis"
    else
        while read line
        do case=$(echo $line | cut -f3 -d" ")
        caseR=$(echo $line | cut -f4 -d" ")
        ctrl=$(echo $line | cut -f1 -d" ")
        ctrlR=$(echo $line | cut -f2 -d" ")
        paste <(cut -f1-5 $ctrl.per_gene_junction_breakdown.txt) <(awk -v n=$ctrlR 'NR==1{for(i=6;i<=NF-1;i++){printf "ctrl_"$i"\t"}{print "ctrl_"$NF}}NR>1{for(i=6;i<=NF-1;i++){if(i==6 || i==8 || i==10){printf $i"\t"}else{printf $i/n"\t"}}{print $i/n}}' $ctrl.per_gene_junction_breakdown.txt) <(awk -v n=$caseR 'NR==1{for(i=6;i<=NF-1;i++){printf "case_"$i"\t"}{print "case_"$NF}}NR>1{for(i=6;i<=NF-1;i++){if(i==6 || i==8 || i==10){printf $i"\t"}else{printf $i/n"\t"}}{print $i/n}}' $case.per_gene_junction_breakdown.txt) | join -t$'\t' -2 1 -1 2 -o 1.1,1.2,1.3,1.4,1.5,2.10,2.11,2.25,2.26,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17 - <(sort -k1,1 $ffDir/$ctrl\_data.normalized.txt | sed 's|total_gene_reads|ctrl_total_gene_reads|;s|FPKM|ctrl_FPKM|;s|total_exon_reads|ctrl_total_exon_reads|;s|total_intron_reads|ctrl_total_intron_reads|;s| \t|\t|g') | join -t$'\t' -2 1 -1 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.10,2.11,2.25,2.26,1.16,1.17,1.18,1.19,1.20,1.21 - <(sort -k1,1 $ffDir/$case\_data.normalized.txt | sed 's|total_gene_reads|case_total_gene_reads|;s|FPKM|case_FPKM|;s|total_exon_reads|case_total_exon_reads|;s|total_intron_reads|case_total_intron_reads|;s| \t|\t|g') > $case\_vs_$ctrl.junction_info_with_multiskipping.txt
	python ~/Desktop/python_scripts/convert_tab_delimited_text_file_to_xlsx.py $case\_vs_$ctrl.junction_info_with_multiskipping.txt $case\_vs_$ctrl.junction_info_with_multiskipping.xlsx
        done < $compFile
fi
