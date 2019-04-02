import re
import sys

checkID = []

with open( sys.argv[1] ) as idFile:
    for i in idFile:
        line = i.strip()
        checkID.append(line)

with open( sys.argv[2] ) as gffFile:
    for i in gffFile:
        line = i.strip()
	line = line.split("\t")
        geneID = line[8].split("=")[2]
	if geneID not in checkID:
            print('\t'.join(map(str,line)))
        else:
            continue
