import sys

checkID = ""
checkGene = ""

with open( sys.argv[1] ) as bedFile:
	for i in bedFile:
		line = i.strip().split("\t")
		if line[20] != ".":
			juncID = line[3]
			gene = line[20].split(":")[1]
			count = line[4]
			#print juncID + "\t" + gene
			if juncID == checkID:
				if checkGene != "" and gene != checkGene:
					print juncID + "\t" + gene + "\t" + checkGene + "\t" + count
			else:
				checkID = juncID
				checkGene = gene
