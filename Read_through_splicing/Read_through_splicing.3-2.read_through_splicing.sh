#!/bin/bash

## This script gets all of the run-on splicing between genes
# input parameters are as follows:
Sample=$1 # first give the sample name
MMR=$2 # mapped reads, in millions
FeatureFolder=$3 # feature file location (file should be labeled $Sample\_data.normalized.txt)
JunctionFolder=$4 # junction file location (file should be labeled $Sample.junctions.bed)
Genome=$5 # genome build
RefFolder=$6 # reference location

# get the overlap of junctions
bedtools intersect -a $JunctionFolder/$Sample.junctions.bed -b $RefFolder/$Genome\_only_exons_refGene_longest.gff3 -wao -split -s > $Sample.intersectJunc.bed

# remove those that appear more or less than twice
# more than twice are isoform issues, not real run-on sequencing/splicing
# less than twice are not in genes
awk 'NR==FNR{junc[$4]++;next}{if(junc[$4]==2){print $0}}' $Sample.intersectJunc.bed $Sample.intersectJunc.bed > $Sample.intersectJunc.clean.bed

# get all of the run-on sequencing/trans splicing
python clean_intersectBed_for_trans_splicing.py $Sample.intersectJunc.clean.bed > $Sample.intersectJunc.ReadThrough.txt

# add in gene names
join -t$'\t' -1 2 -2 1 -o 1.1,1.2,2.2,1.3,1.4 <(sort -k2,2 $Sample.intersectJunc.ReadThrough.txt) <(sort -k1,1 $RefFolder/$Genome\_gene_name_ID) | sort -k4,4 | join -t$'\t' -1 4 -2 1 -o 1.1,1.2,1.3,1.4,2.2,1.5 - <(sort -k1,1 $RefFolder/$Genome\_gene_name_ID) > $Sample.intersectJunc.gene_names.txt

# clean this up
awk -v m=$MMR '{print $3":"$5"\t"$6/m}' $Sample.intersectJunc.gene_names.txt | awk '{sum[$1]+=$2}END{for(i in sum){print i"\t"sum[i]}}' | sed 's|:|\t|g' > $Sample.intersectJunc.gene_names.clean.txt

# add in FPKM and normalize the junction reads to RPM
join -t$'\t' -1 1 -2 2 -o 1.1,2.1,2.3,2.10,2.11,2.38,1.2,1.3 <(sort -k1,1 $Sample.intersectJunc.gene_names.clean.txt) <(tail -n+2 $FeatureFolder/$Sample\_data.normalized.txt | sort -k2,2) | sort -k7,7 | join -t$'\t' -1 7 -2 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.1,2.3,2.10,2.11,2.38,1.8 - <(tail -n+2 $FeatureFolder/$Sample\_data.normalized.txt | sort -k2,2) > $Sample.ReadThroughJunctions.all.txt

python remove_overlaps_get_intergene_distance.py $RefFolder/$Genome\_mrna_longest.gff3 $Sample.ReadThroughJunctions.all.txt $Sample.ReadThroughJunctions.noOverlap.txt

grep -v -e "SCARNA" -e "SNOR" -e "MIR" -e "POLR" -e "KIAA" -e "LINC" -e "RN7S" -e "FAM" -e "RRP" -e "HIST" -e "orf" $Sample.ReadThroughJunctions.noOverlap.txt | awk '{if($1!~"-"){print $0}}' | awk '{if($6!~"-"){print $0}}' | awk '{if($1!~/^LOC/){print $0}}' | awk '{if($1!~/^RNA/){print $0}}' > $Sample.ReadThroughJunctions.filtered.txt


# remove the temporary files
rm $Sample.intersectJunc.bed
rm $Sample.intersectJunc.clean.bed
rm $Sample.intersectJunc.ReadThrough.txt
rm $Sample.intersectJunc.gene_names.txt
rm $Sample.intersectJunc.gene_names.clean.txt
