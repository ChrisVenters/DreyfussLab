import re
import sys

checkID = {}

with open( sys.argv[1] ) as idFile:
    for i in idFile:
        line = i.strip().split("\t")
        inputID = line[0]
        if int(line[2]) >= -500:
		continue
	elif int(line[2]) >= -5000:
        	intDistance = int(line[2])
		checkID[inputID] = intDistance
        else:
		intDistance = -5000
		checkID[inputID] = intDistance

with open( sys.argv[2] ) as gffFile:
    for i in gffFile:
        line = i.strip()
        line = line.split("\t")
        geneID = line[8].split("=")[1].split(";")[0]
        geneName = line[8].split("=")[2]
        idField = "ID=" + geneID + "_upstream;Parent=" + geneName + "_upstream"
        if geneID in checkID:
		if line[6] == "+":
			geneEnd = int(line[3])
			geneStart = geneEnd + int(checkID[geneID])
			strand = "-"
		else:
			geneStart = int(line[4])
			geneEnd = geneStart - int(checkID[geneID])
			strand = "+"
		print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line[0], line[1], line[2], geneStart, geneEnd, line[5], strand, line[7], idField)
        else:
            continue
