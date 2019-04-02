## use last CDS to replace last exon when the last exon has a CDS region

import re
import sys

dict_maxExon = {}
f1 = open(sys.argv[1])
for i in f1:
	l = i.strip().split('\t')
	exon_num = l[4]
	dict_maxExon[l[1]+'\t'+exon_num] = ''
f1.close()

## last exon coors
dict_exon_coor = {}
f2 = open(sys.argv[2])
for j in f2:
	l2 = j.strip().split('\t')
	refseqid = l2[8].strip().split(';')[0].split(':')[1]
	exon_num2 = l2[8].strip().split(';')[0].split(':')[2]
	if re.search('exon',j):
		if refseqid+'\t'+exon_num2 in dict_maxExon:
			start = l2[3] 
			end = l2[4]
			strand = l2[6]
			dict_exon_coor[refseqid] = l2[0]+"\t"+refseqid+'\t'+start+'\t'+end+"\t"+strand+'\t'+exon_num2

f2.close()

## last exon CDS reads
dict_lecds = {}
f3 = open(sys.argv[2])
for k in f3:
	if re.search('CDS',k):
		l3 = k.strip().split('\t')
		refseqid2 = l3[8].strip().split(';')[0].split(':')[1]
		if refseqid2 in dict_exon_coor:
			if l3[6] == "+":
			## last exon CDS and its corresponding last exon (id+number)
				if int(l3[3]) >= int(dict_exon_coor[refseqid2].split('\t')[2]):
					lecds = refseqid2 + '\t' + dict_exon_coor[refseqid2].split('\t')[5]
					dict_lecds[lecds] = l3[9]
			else:
				if int(l3[4]) <= int(dict_exon_coor[refseqid2].split('\t')[3]):
					lecds = refseqid2 + '\t' + dict_exon_coor[refseqid2].split('\t')[5]
					dict_lecds[lecds] = l3[9]

f3.close()

## output format
for each in sys.stdin:
	line = each.strip().split('\t')
	key =  line[0]+'\t'+ str(int(line[1])+1)
	if key in dict_lecds:
		print '\t'.join(line[0:5]) + '\t' + dict_lecds[key] + '\t'+ line[6] + '\t' + line[7]

	else:
		# exclude the intron which nextExon are pure 3'UTR
		if key in dict_maxExon:
			pass
		else:
			print each.strip()
