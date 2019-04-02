#!/bin/bash

SAMPLE=$1 # sample name
GENOME=$2 # genome build
WHOLEGENEPCPA=$3 # whole gene PCPA file
EXPRESSIONCHANGE=$4 # directory of expression change files
INTRONICPCPA=$5 # intron based PCPA file
RANDOMSAMPLING=$6 # text file of distances between randomly sampled genes
# note, the random sampling I did separately for each organism and the files are: ~/Desktop/pipelines/gene_neighborhoods/

# run PCPA calculations for whole gene and intron

# run up/down regulation

# run 3'UTR shortening

# do the neighborhood analysis
cut -f2 $WHOLEGENEPCPA | tail -n+2 > $SAMPLE.pcpaed_genes.IDs
cut -f1 $EXPRESSIONCHANGE/$SAMPLE.no_exon1_up.txt | tail -n+2 > $SAMPLE.no_exon1_upregulated.IDs
cut -f1 $EXPRESSIONCHANGE/$SAMPLE.no_exon1_down.txt | tail -n+2 > $SAMPLE.no_exon1_downregulated.IDs
cut -f1 $EXPRESSIONCHANGE/$SAMPLE.all_exons_down.txt | tail -n+2 > $SAMPLE.downregulated.IDs


for file in *.IDs; do n=$(echo $file | sed 's|IDs|gff3|'); while read line; do grep -w $line ~/Desktop/pipelines/organism_lib/$GENOME.refGene_gene.gff3; done < $file | sort -k1,1 -k4,4n > $n; done

closestBed -io -d -a $SAMPLE.no_exon1_upregulated.gff3 -b $SAMPLE.pcpaed_genes.gff3 > $SAMPLE.up_to_pcpa.bed
closestBed -io -d -a $SAMPLE.no_exon1_upregulated.gff3 -b $SAMPLE.no_exon1_downregulated.gff3 > $SAMPLE.up_to_down_no_exon1.bed
closestBed -io -d -a $SAMPLE.no_exon1_upregulated.gff3 -b $SAMPLE.downregulated.gff3 > $SAMPLE.up_to_down.bed
closestBed -io -d -a $SAMPLE.no_exon1_upregulated.gff3 -b $SAMPLE.no_exon1_upregulated.gff3 > $SAMPLE.up_to_up.bed

# clean up the distances for plotting
for file in $SAMPLE.*.bed; do n=$(echo $file | sed 's|.bed|.txt|'); awk '{print $NF}' $file > $n; done

if [ ! -e "intermediate" ]
then
	mkdir intermediate
fi

# with intronic PCPA
if [ ! -e "with_intronic_pcpa" ]
then
	mkdir with_intronic_pcpa
fi
cd with_intronic_pcpa/
cut -f2 $INTRONICPCPA | tail -n+2 | awk '!a[$1]++' > $SAMPLE.intronic_pcpa.IDs
while read line; do grep -w $line ~/Desktop/pipelines/organism_lib/$GENOME.refGene_gene.gff3; done < $SAMPLE.intronic_pcpa.IDs | sort -k1,1 -k4,4n > $SAMPLE.intronic_pcpa.gff3
cat $SAMPLE.intronic_pcpa.IDs ../$SAMPLE.pcpaed_genes.IDs | awk '!a[$NF]++' > $SAMPLE.combined_pcpaed_genes.IDs
cat $SAMPLE.intronic_pcpa.gff3 ../$SAMPLE.pcpaed_genes.gff3 | awk '!a[$NF]++' | sort -k1,1 -k4,4n > $SAMPLE.combined_pcpaed_genes.gff3
closestBed -io -d -a ../$SAMPLE.no_exon1_upregulated.gff3 -b $SAMPLE.combined_pcpaed_genes.gff3 > $SAMPLE.up_to_pcpa.bed
closestBed -io -d -a $SAMPLE.combined_pcpaed_genes.gff3 -b ../$SAMPLE.no_exon1_upregulated.gff3 > $SAMPLE.pcpa_to_up.bed

comm -13 <(sort $SAMPLE.combined_pcpaed_genes.IDs) <(sort ../$SAMPLE.downregulated.IDs) > $SAMPLE.down_not_PCPAed.IDs
while read line; do grep -w $line ~/Desktop/pipelines/organism_lib/$GENOME.refGene_gene.gff3; done < $SAMPLE.down_not_PCPAed.IDs | sort -k1,1 -k4,4n > $SAMPLE.down_not_PCPAed.gff3
closestBed -io -d -a ../$SAMPLE.no_exon1_upregulated.gff3 -b $SAMPLE.down_not_PCPAed.gff3 > $SAMPLE.up_to_only_down.bed
closestBed -io -d -a $SAMPLE.down_not_PCPAed.gff3 -b ../$SAMPLE.no_exon1_upregulated.gff3 > $SAMPLE.only_down_to_up.bed

# clean up the distances for plotting
for file in $SAMPLE.*.bed; do n=$(echo $file | sed 's|.bed|.txt|'); awk '{print $NF}' $file > $n; done
# prep data to plot a histogram in R for the frequency of distance between genes

if [ ! -e "intermediate" ]
then
	mkdir intermediate
fi

mv $SAMPLE.*.bed intermediate/
mv $SAMPLE.*.gff3 intermediate/
mv $SAMPLE.*.IDs intermediate/

# for boxplot
cat <(echo -e "Group\tDistance") <(awk '{print "Random\t"$1}' $RANDOMSAMPLING) <(awk '{print "Up_to_PCPA\t"$1}' $SAMPLE.up_to_pcpa.txt) <(awk '{print "Up_to_Down\t"$1}' $SAMPLE.up_to_only_down.txt) <(awk '{print "Up_to_Up\t"$1}' ../$SAMPLE.up_to_up.txt) > ../$SAMPLE.group_breakdown_distance.boxplot.txt

cd ../

mv $SAMPLE.*.bed intermediate/
mv $SAMPLE.*.gff3 intermediate/
mv $SAMPLE.*.IDs intermediate/

Rscript ~/Desktop/pipelines/gene_neighborhoods/makePlots_noUTRShortening.R $SAMPLE.group_breakdown_distance.boxplot.txt $SAMPLE
