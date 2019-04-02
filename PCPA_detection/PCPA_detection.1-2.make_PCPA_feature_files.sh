#!/bin/bash

# creates the files needed to calculate PCPA

sample=$1 # ex. control_serum
refPath=$2 # ex. /home/venters/pipelines/organism_lib/
genome=$3 # ex. rn6
bamFile=$4 # the cleaned bam file
libPath=$5 # intron/exon features analyzed for PCPA 
featureFile=$6 # this is the raw reads feature file I make, ex. /home/venters/projects/qiang_rat_cells/feature_files/control_serum_data.raw_reads.txt
rpkmFile=$7 # this is exon count normalized GFold file, ex. /home/venters/projects/qiang_rat_cells/GFold/total_mapped_reads_rpkm/control_serum.read_count

mapDepth=$(samtools idxstats $bamFile | awk '{sum+=$3}END{print sum}')

### first get the coverage needed
coverageBed -abam $bamFile -b $refPath/$genome.longest_isoform_splitIntron_Exon.gff3 -split > $sample.splitIntron_Exon.bed
coverageBed -abam $bamFile -b $refPath/$genome.longest_isoform_splitIntron1stlastExon.gff3 -split > $sample.splitIntron1stlastExon.bed
coverageBed -abam $bamFile -b $refPath/$genome.longest_isoform_utrCDSExon.gff3 -split > $sample.longest_isoform_utrCDSExon.bed

### Now the actual file making

# format each intron, add intron order, and reads for intron1 quarters, previous exon, next exon 
awk -v n=$sample '{if($3=="intron_quarter1"){print $0 >> n".intronQ1"}else if($3=="intron_quarter4"){print $0 >> n".intronQ4"}else if($3=="exon"){print $0 > n".exon"}}' $sample.splitIntron1stlastExon.bed
### THIS IS THE SAME AS CHAO'S
python $libPath/format.py $sample.intronQ1 $sample.intronQ4 $sample.exon | sort -k1,1 > $sample.format1

# step2, replace last exon with last cds, if its lastexon has no CDS(pure 3UTR,6000 genes), keep it at it is
cat $sample.format1 | python $libPath/lastExon2CDS.py $refPath/$genome.longest_isoform.txt $sample.longest_isoform_utrCDSExon.bed | sort -k1,1 -k2,2n | cat <(echo -e "refseqID\tquarter1\tquarter4\tfirstExon\tlastExon") - > $sample.format2

# step3, add features, gene name, length etc, add read count, rpkm etc and exon-exon junction reads
cat $sample.format2 | python $libPath/addfeature.py $refPath/$genome.longest_isoform.txt $rpkmFile $featureFile $mapDepth | sort -k1,1 | cat <(echo -e "gene_name\trefseqID\tgene_length\ttranscript_length\ttotal_exon_number\trpkm\tquarter1\tnorm_quarter1\tquarter4\tnorm_quarter4\tfirstExon\tnorm_firstExon\tlastExon\tnorm_lastExon\texon_exon_junc\tnorm_exon_exon_junc\ttotal_exon\tnorm_total_exon\ttotal_intron\tnorm_total_intron") - > $sample.feature
