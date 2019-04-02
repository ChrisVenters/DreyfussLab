## add more features for PCPA calculations

import sys
import re

## use the features in longest_isoform.txt,gene_name\trefseqID\tgene_length\ttranscript_length\ttotal_exon_number 
dict1 = {}
f1 = open(sys.argv[1])
for i in f1:
	l = i.strip().split('\t')
	refseqid = l[1]
	dict1[refseqid] = i.strip()
f1.close()

## use the features in longest_isoform_intronExon.gff3, intron_order and intron_length
dict2 = {}
f2 = open(sys.argv[2])
for j in f2:
	if re.search('intron',j):
		l2 = j.strip().split('\t')
		refseqid2 = l2[8].strip().split(';') [1].split('=')[1]
		intron_order = l2[8].strip().split(';') [0].split(':')[2] 
		intron_len = int(l2[4]) - int(l2[3]) + 1
		dict2[refseqid2+'\t'+intron_order] = intron_len
f2.close()

## Gfold reads count file, only use rpkm
dict3 = {}
f3 = open(sys.argv[3])
for k in f3:
	l3 = k.strip().split('\t')
	geneName = l3[0]
	dict3[geneName] = l3[4]
f3.close()


## output format
for each in sys.stdin:
	if not re.match('refseqID',each):
		line = each.strip().split('\t')	
		intron_order = line[1]
		total_exon_reads = line[6]
		total_intron_reads = line[7]
		## read counts normalize to total mapped reads, get reads per million
		norm_q1 = round(float(line[2])/int(sys.argv[4])*1000000,4) 
		norm_q4 = round(float(line[3])/int(sys.argv[4])*1000000,4)
		norm_pE = round(float(line[4])/int(sys.argv[4])*1000000,4)
		norm_nE = round(float(line[5])/int(sys.argv[4])*1000000,4)
		norm_intron =  round(float(total_intron_reads)/int(sys.argv[4])*1000000,4)
		if line[0] in dict1:
			prefix1 = dict1[line[0]]
			gene_name = prefix1.split('\t')[0]
		# add intron length
		if line[0]+'\t'+line[1] in dict2:
			prefix2 = dict2[line[0]+'\t'+line[1]]
		if gene_name in dict3:	
			prefix3 = dict3[gene_name]
		
		print prefix1 +'\t'+intron_order+'\t'+str(prefix2)+'\t'+total_exon_reads+'\t'+prefix3+'\t'+total_intron_reads+'\t'+str(norm_intron)+'\t'+line[2]+'\t'+ str(norm_q1)+'\t'+line[3]+'\t'+str(norm_q4)+'\t'+line[4]+'\t'+str(norm_pE)+'\t'+line[5]+'\t'+str(norm_nE)





