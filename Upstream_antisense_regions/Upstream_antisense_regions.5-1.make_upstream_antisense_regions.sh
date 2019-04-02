#!/bin/bash

GENOME=$1

# get a working file of all transcripts
awk '$3=="transcript" || $3=="mRNA"' ~/Desktop/pipelines/organism_lib/$GENOME.refGene_full.gff3 > $GENOME.refGene_genes.gff3

# get rid of miRNAs, snoRNAs, scaRNAs, AS sequences
awk '{if($3=="transcript"){print $0}}' $GENOME.refGene_genes.gff3 | cut -f3 -d"=" | awk '/^MIR/ || /^SNO/ || /^SCARNA/ || /LOC/ || /HNRNP/ || /^RNU/ || /^SNAR/ || /RNV/ || /^LINC/ || /^FAM/ || /^ANKRD/ || /^ZNF/ || $1~"-AS" || $1~"-IT" || /orf/ || /^Mir/ || /^Sno/ || /^Scarna/ || /^Rnu/ || /^Linc/ || /Rik$/ || /^Fam/ || /^mir/ || /^sno/ || /^scaRNA/ || /^fam/ || /:/ || /-/ || /^Gm/ || /CR[0-9]/ || /-rs/ || /-ps/ || /rRNA/ || /snRNA/ || /snmRNA/ || /SrRNA/ || /7sk/ || /7sl/ || /rnu/' > $GENOME.refGene_junk_transcripts.txt

python remove_junk_transcript.py $GENOME.refGene_junk_transcripts.txt $GENOME.refGene_genes.gff3 > $GENOME.refGene_no_junk_transcripts.gff3

# get out the longest isoform from the gene list that doesn't include ncRNAs
cut -f1,2 ~/Desktop/pipelines/organism_lib/$GENOME.longest_isoform.txt > $GENOME.refGene_longest_isoform.txt
cut -f2 $GENOME.refGene_longest_isoform.txt > $GENOME.refGene_longest_isoform.ID

while read line; do grep -w $line $GENOME.refGene_no_junk_transcripts.gff3; done < $GENOME.refGene_longest_isoform.ID > $GENOME.refGene_clean_genes.gff3

sort -k1,1 -k4,4n $GENOME.refGene_clean_genes.gff3 > tmp.gff3; mv tmp.gff3 $GENOME.refGene_clean_genes.gff3
sort -k1,1 -k4,4n $GENOME.refGene_genes.gff3 > tmp.gff3; mv tmp.gff3 $GENOME.refGene_genes.gff3

# dont double dip in genes that are upstream antisense from one another
# first check for genes <=5kb from one another
# do this 2 ways: 1) check against our clean gene list
closestBed -S -D a -id -a $GENOME.refGene_clean_genes.gff3 -b $GENOME.refGene_clean_genes.gff3 > $GENOME.check_upstream_genes.txt
# 2) also do it without excluding ncRNAs
closestBed -S -D a -id -a $GENOME.refGene_clean_genes.gff3 -b $GENOME.refGene_genes.gff3 > $GENOME.check_upstream_genes_with_ncRNAs.txt

# there appear to be a couple of thousand cases where two different strand genes overlap on their 3' ends (as in, they're NOT upstream antisense) and BEDTools thinks overlap means genes are both upstream and downstream from one another
# so I need to remove these cases myself from the overlap and re-run them to check for upstream antisense genes
# 1)
awk '{if($7=="+" && $14>$4 && $18!="."){print $0}}' $GENOME.check_upstream_genes.txt | awk '!a[$9]++' | cut -f1-9 > $GENOME.incorrect_overlap_genes.gff3
awk '{if($7=="-" && $5>$13 && $18!="."){print $0}}' $GENOME.check_upstream_genes.txt | awk '!a[$9]++' | cut -f1-9 >> $GENOME.incorrect_overlap_genes.gff3
sort -k1,1 -k4,4n $GENOME.incorrect_overlap_genes.gff3 > tmp.gff3; mv tmp.gff3 $GENOME.incorrect_overlap_genes.gff3
closestBed -S -D a -id -io -a $GENOME.incorrect_overlap_genes.gff3 -b $GENOME.refGene_clean_genes.gff3 > $GENOME.check_upstream_genes_fixed_incorrect_overlap.txt
# remove the incorrect genes from the base file
awk '{if($7=="+" && $14>$4){next}else{print $0}}' $GENOME.check_upstream_genes.txt | awk '{if($7=="-" && $5>$13){next}else{print $0}}' > $GENOME.check_upstream_genes.clean.txt

# 2)
awk '{if($7=="+" && $14>$4 && $18!="."){print $0}}' $GENOME.check_upstream_genes_with_ncRNAs.txt | awk '!a[$9]++' | cut -f1-9 > $GENOME.incorrect_overlap_genes_with_ncRNAs.gff3
awk '{if($7=="-" && $5>$13 && $18!="."){print $0}}' $GENOME.check_upstream_genes_with_ncRNAs.txt | awk '!a[$9]++' | cut -f1-9 >> $GENOME.incorrect_overlap_genes_with_ncRNAs.gff3
sort -k1,1 -k4,4n $GENOME.incorrect_overlap_genes_with_ncRNAs.gff3 > tmp.gff3; mv tmp.gff3 $GENOME.incorrect_overlap_genes_with_ncRNAs.gff3
closestBed -S -D a -id -io -a $GENOME.incorrect_overlap_genes_with_ncRNAs.gff3 -b $GENOME.refGene_clean_genes.gff3 > $GENOME.check_upstream_genes_fixed_incorrect_overlap_with_ncRNAs.txt
# remove the incorrect genes from the base file
awk '{if($7=="+" && $14>$4){next}else{print $0}}' $GENOME.check_upstream_genes_with_ncRNAs.txt | awk '{if($7=="-" && $5>$13){next}else{print $0}}' > $GENOME.check_upstream_genes_with_ncRNAs.clean.txt


if [ ! -e "with_ncrnas" ]
then
	mkdir with_ncrnas
fi

mv *_with_ncRNAs* with_ncrnas/

# take only one occurence of each gene, but the first occurence should be the one with the smallest intergenic distance
cut -f9,18,19 $GENOME.check_upstream_genes.clean.txt | awk '{if($2=="."){print $1"\t-10000"}else{print $1"\t"$3}}' | sort -rnk2,2 | awk '!a[$1]++'| sed 's|ID=||;s|;Parent=|\t|' > $GENOME.genes_with_upstream_distances.txt
# then create upstream antisense gene intervals 5kb or less that don't overlap with known genes
python make_upstream_regions.py $GENOME.genes_with_upstream_distances.txt $GENOME.refGene_genes.gff3 > $GENOME.upstream_regions.gff3
sort -k1,1 -k4,4n $GENOME.upstream_regions.gff3 > tmp.gff3; mv tmp.gff3 $GENOME.upstream_regions.gff3
awk '{if($4<0 || $5 <0){next}else{print $0}}' $GENOME.upstream_regions.gff3 > tmp.gff3; mv tmp.gff3 $GENOME.upstream_regions.gff3
