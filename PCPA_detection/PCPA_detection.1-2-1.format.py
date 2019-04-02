## fortmat the feature reads count for statistic tests

import sys
import re

dict = {}
f1 = open(sys.argv[1])
for i in f1:
	if re.search("intron_quarter1",i):
		l = i.strip().split('\t')
		anno = l[8].strip().split(';')[0] 
		ID = anno.split(":")[1]
		intron_num = anno.split(":")[2] 
		dict[ID+'\t'+intron_num] = [l[9]]
f1.close()

f2 = open(sys.argv[2])
for j in f2:	
	if re.search("intron_quarter4",j):
		l2 = j.strip().split('\t')
		anno2 = l2[8].strip().split(';')[0] 
		ID2 = anno2.split(":")[1]
		intron_num2 = anno2.split(":")[2]
		key2 = ID2 + '\t' + intron_num2
		if key2 in dict:		
			dict[key2].append(l2[9])
	
f2.close()


f3 = open(sys.argv[3])
for k in f3:
	l3 = k.strip().split('\t')
	anno3 = l3[8].strip().split(';')[0] 
	ID3 = anno3.split(":")[1]
	if re.search("exon",k):
		exon_num = int(anno3.split(":")[2])
		# previous exon,same order with this intron
		key3 = ID3 + '\t' + str(exon_num)
		if key3 in dict:
			dict[key3].append(l3[9])
f3.close()

f4 = open(sys.argv[3])
for m in f4:
	l4 = m.strip().split('\t')
	anno4 = l4[8].strip().split(';')[0] 
	ID4 = anno4.split(":")[1]
	if re.search("exon",m):
		exon_num2 = int(anno4.split(":")[2])
		# next exon
		key4 = ID4 + '\t' + str(exon_num2 - 1)
		if key4 in dict:
			dict[key4].append(l4[9])
f4.close()





#print "refseqID\tintron_num\tquater1\tquater4\tnextExon"
for each in dict:
	print each + '\t'+ '\t'.join(dict[each])



