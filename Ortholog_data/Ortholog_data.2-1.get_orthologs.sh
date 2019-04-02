# download Ensembl ortholog data from Biomart
# human_ensembl_gene_ID, specis_ensembl_transcript_ID, gene_start, gene_end
# ex. ensembl.mouse_human_orthologues.txt
# also download Ensembl IDs to gene name
# ex. human_ensembl_to_gene_name.txt

# download ensgene from UCSC table browser (ensGene format)
# ex. ensgene.mouse.txt

## First get the genes affected by U1 AMO from human
# extract relevant (ex. PCPAed) genes from Ensembl data
join -t$'\t' -1 3 -2 1 <(tail -n+2 human_ensembl_to_gene_name.txt | sort -k3,3) <(cut -f1 U1_AMO_8hr_labeled.pcpa.filter.txt | tail -n+2 | sort) | awk '!a[$2]++' > U1_AMO_8hr_labeled.pcpa.filter.ensembl_id

# get the gene lengths
while read line; do n=$(echo $line | cut -f1 -d" "); m=$(echo $line | cut -f2 -d" "); grep -w $m ensgene.human.txt | awk '{print $6-$5}' | sort -rn | head -n1 | paste <(echo $n) <(echo $m) -; done < U1_AMO_8hr_labeled.pcpa.filter.combined.ensembl_id | sort -rnk3,3 | awk '!a[$1]++' | sort -k1,1 > U1_AMO_8hr_labeled.pcpa.filter.combined.gene_length


# do the same for up-regulated
join -t$'\t' -1 3 -2 1 <(tail -n+2 human_ensembl_to_gene_name.txt | sort -k3,3) <(sort non-PCPAedgene_exons_plus_intronless_up.id) | awk '!a[$2]++' > non-PCPAedgene_exons_plus_intronless_up.ensembl_id

while read line; do n=$(echo $line | cut -f1 -d" "); m=$(echo $line | cut -f2 -d" "); grep -w $m ensgene.human.txt | awk '{print $6-$5}' | sort -rn | head -n1 | paste <(echo $n) <(echo $m) -; done < non-PCPAedgene_exons_plus_intronless_up.combined.ensembl_id | sort -rnk3,3 | awk '!a[$1]++' | sort -k1,1 > non-PCPAedgene_exons_plus_intronless_up.combined.gene_length

## Now get orthologs from other organisms
# first get all orthologs
while read line; do n=$(echo $line | cut -f1 -d" "); join -t$'\t' <(tail -n+2 human_ensembl_to_gene_name.combined.txt | sort -k1,1 | cut -f1-3) <(tail -n+2 ensembl.$n\_human_orthologues.txt | sort -k1,1 | cut -f1,2) > all_human_overlap.$n.txt; done < species.txt
# unique results only
while read line; do n=$(echo $line | cut -f1 -d" "); m=$(echo $line | cut -f2 -d" "); join -t$'\t' -1 4 -2 2 <(sort -k4,4 all_human_overlap.$n.txt) <(tail -n+2 $m\_ensembl.txt | sort -k2,2) | awk '!a[$NF]++' > $m.human_orthologs.uniq_gene_names.txt; done < species.txt

# overlap with PCPAed or upregulated genes
while read line; do n=$(echo $line | cut -f1 -d" "); m=$(echo $line | cut -f2 -d" "); join -t$'\t' <(cut -f4,6 $n.human_orthologs.uniq_gene_names.txt | sort -k1,1) <(sort -k1,1 PCPA_data/U1_AMO_8hr_labeled.pcpa.filter.txt | cut -f1) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(tail -n+2 PCPA_data/$m.pcpa_filter.txt | sort -k1,1) | cat <(echo -e "gene_name\thuman_gene_name\tgene_ID\tgene_length\ttranscript_length\ttotal_exon_number\tctrl_lastExon_VS_u1amo_lastExon\tpoissonP\tfisherP\tpoissonQ\tfisherQ") - > ortholog_data/$m.pcpa_filter.human_orthologs.txt; join -t$'\t' <(cut -f4,6 $n.human_orthologs.uniq_gene_names.txt | sort -k1,1) <(sort -k1,1 PCPA_data/U1_AMO_8hr_labeled.no_pcpa_exons_up.txt | cut -f1) | sort -k2,2 | join -t$'\t' -1 2 -2 2 - <(tail -n+2 PCPA_data/$m.no_pcpa_exons_up.txt | sort -k2,2) | cat <(echo -e "gene_name\thuman_gene_name\tgene_ID\tgene_length\tmRNA_length\tctrl_no_exon1_VS_case_no_exon1_RPM\tctrl_all_exons_VS_case_all_exons_RPM\tcase_no_exon1_VS_ctrl_no_exon1_RPM\tcase_all_exons_VS_ctrl_all_exons_RPM\texons_up.q\texons_down.q") - > ortholog_data/$m.no_pcpa_exons_up.human_orthologs.txt; done < genomes.txt

ls ortholog_data/*pcpa_filter*.txt | sed 's|ortholog_data/||;s|.pcpa_filter.human_orthologs.txt||' > organism_samples.txt

while read line; do cat <(awk 'NR==1{print $3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\tortholog_status"}NR>1{print $3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\toverlap"}' ortholog_data/$line.no_pcpa_exons_up.human_orthologs.txt) <(awk 'NR==FNR{a[$3];next} !($1 in a)' ortholog_data/$line.no_pcpa_exons_up.human_orthologs.txt PCPA_data/$line.no_pcpa_exons_up.txt | awk '{print $0"\tno_overlap"}') > for_graphs/$line.no_pcpa_exons_up.ortho_status.txt; cat <(awk 'NR==1{print $3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\tortholog_status"}NR>1{print $3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\toverlap"}' ortholog_data/$line.pcpa_filter.human_orthologs.txt) <(awk 'NR==FNR{a[$3];next} !($2 in a)' ortholog_data/$line.pcpa_filter.human_orthologs.txt PCPA_data/$line.pcpa_filter.txt | awk '{print $0"\tno_overlap"}') > for_graphs/$line.pcpa_filter.ortho_status.txt; done < organism_samples.txt
