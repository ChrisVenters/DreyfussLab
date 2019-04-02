import sys
# this script will remove overlapping genes and print out the distance between the genes of those that remain

geneLoc = {}
outputLines = []

# load the gene starts and ends into the array 
# use gff3 mRNA file (ex. ~/Desktop/pipelines/organism_lib/hg38_mrna_longest.gff3)
with open( sys.argv[1] ) as gff3File:
	for i in gff3File:
		line = i.strip().split("\t")
		chrom = line[0]
		strand = line[6]
		# need to change the start and end based on strand
		start = int(line[3])
		end = int(line[4])
		geneID = line[8].split(";")[0].split("=")[1]
		# make a dict in a dict for the gene info
		geneLoc[geneID] = {}
		geneLoc[geneID]["chrom"] = chrom
		geneLoc[geneID]["start"] = start
		geneLoc[geneID]["end"] = end
		geneLoc[geneID]["strand"] = strand
# load the file from the current script
# here we'll get the gene IDs for the trans-splicing and check against the previous array
with open( sys.argv[2] ) as geneFile:
	for i in geneFile:
		line = i.strip().split("\t")
		geneOne = line[1]
		geneTwo = line[7]
		# first find out if the genes are overlapping, here are the possibilities
		# 1) gene one starts inside gene two
		if geneLoc[geneOne]["start"] >= geneLoc[geneTwo]["start"] and geneLoc[geneOne]["start"] <= geneLoc[geneTwo]["end"]:
			continue
		# 2) gene two starts inside gene one
		elif geneLoc[geneTwo]["start"] >= geneLoc[geneOne]["start"] and geneLoc[geneTwo]["start"] <= geneLoc[geneOne]["end"]:
			continue
		# 3) gene one ends inside gene two
		elif geneLoc[geneOne]["end"] <= geneLoc[geneTwo]["end"] and geneLoc[geneOne]["end"] >= geneLoc[geneTwo]["start"]:
			continue
		# 4) gene two ends inside gene one
		elif geneLoc[geneTwo]["end"] <= geneLoc[geneOne]["end"] and geneLoc[geneTwo]["end"] <= geneLoc[geneOne]["end"]:
			continue
		# if the genes don't overlap, then we compute distance between
		else:
			# we need to determine which gene is the upstream gene
			if geneLoc[geneOne]["start"] < geneLoc[geneTwo]["start"]:
				# then compute distance from 3' end of upstream gene to 5' end of downstream
				intergeneDistance = str(geneLoc[geneTwo]["start"] - geneLoc[geneOne]["end"])
			else:
				intergeneDistance = str(geneLoc[geneOne]["start"] - geneLoc[geneTwo]["end"])
			line.append(intergeneDistance)
			outputLines.append(line)

with open( sys.argv[3], 'w' ) as outFile:
	outFile.write("Gene1\tGene1ID\tLength\tGeneReads\tRPKM\tJunctionReads\tGene2\tGene2ID\tLength\tGeneReads\tRPKM\tJunctionReads\tReadThroughSplicing\tIntergeneDistance\n")
	for x in outputLines:
		outLine = "\t".join(x)
		outFile.write( "%s\n" % ( outLine ) ) 
