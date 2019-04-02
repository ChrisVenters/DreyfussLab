# initial fasta files downloaded from biomart on 1-19-16 from GRCh38.p5
# filter for snRNA or miscRNA sequences
# attributes: Gene Name, Unspliced sequence (transcript), Ensembl Gene ID and Transcript ID

# clean out the 7SK and 7SL from the miscRNA
awk 'BEGIN{RS=">"}{if($0~"7S"){print RS$0}}' miscRNA.biomart.fa | awk '{if($1==""){next}else{print $0}}' | cut -f1,2 -d"|" > miscRNA.biomart.clean.fa

# clean out the RNUX from the snRNA
# then I can make a gtf file with an entry for each snRNA that incorporates all the entries for that snRNA
awk 'BEGIN{RS=">"}{if($0~"U"){print RS$0}}' snRNA.biomart.fa | awk '{if($1==""){next}else{print $0}}' | tr '[:lower:]' '[:upper:]' | cut -f1,2 -d"|" > snRNA.biomart.clean.fa 

# combine the snRNA and miscRNA
cat snRNA.biomart.clean.fa miscRNA.biomart.clean.fa > snRNA.fa

# making the gtf file for use with the metagene
awk 'BEGIN{RS=">";OFS="\t"}{split($1,k,"|")}{j=0;for(i=2;i<=NF;i++){j+=length($i)}{print $1,"snRNA\texon\t1",j,".\t+\t.\tgene_id \""k[1]":"k[2]"\"; transcript_id \""k[2]"\""}}' snRNA.fa | tail -n+2 > snRNA.gtf

## clean up the pseudogenes and link to the ensembl gene ID
# this makes a file for each snRNA (U1, U2, etc.) that contains all the ensembl gene IDs for that snRNA
cut -f1 snRNA.gtf | sed 's/-.*|ENSG/|ENSG/g;s/ATAC.*P/ATAC/g;s|RN||g;s/U5A\|U5B\|U5D\|U5E\|U5F/U5/g;s|VU1|U1|g;s|U6V|U6|g;s/P.*|ENSG/|ENSG/g;s/L.*|ENSG/L|ENSG/g' | awk 'BEGIN{FS="|"}{a[$2]=$1}END{for(i in a){print i >> a[i]".ensembl_id.txt"}}'
 
# make a list of each snRNA
cut -f1 -d"|" snRNA.gtf | sed 's/-.*|ENSG/|ENSG/g;s/ATAC.*P/ATAC/g;s|RN||g;s/U5A\|U5B\|U5D\|U5E\|U5F/U5/g;s|VU1|U1|g;s|U6V|U6|g;s/P.*|ENSG/|ENSG/g;s/L.*|ENSG/L|ENSG/g' | awk '!a[$1]++' > snRNA_list.txt
# separate each snRNA out into it's own GTF for the metagene
while read line; do while read ens; do grep -w $ens snRNA.gtf; done < $line.ensembl_id.txt > $line\_snRNA.gtf; done < snRNA_list.txt

# make a Bowtie2 index
bowtie2-build snRNA.fa snRNA
