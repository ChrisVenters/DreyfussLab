import sys

starts = {}
ends = {}
outJunc = []

with open( sys.argv[1] ) as gff3File:
	for i in gff3File:
		line = i.strip().split("\t")
		chrom = line[0]
		start = line[3]
		end = line[4]
		gene = line[8].split("=")[2]
		cStart = chrom + ":" + start
		cEnd = chrom + ":" + end
		starts[cStart] = gene
		ends[cEnd] = gene

with open( sys.argv[2] ) as juncFile:
	for i in juncFile:
		line = i.strip().split("\t")
		chrom = line[0]
		start = line[1]
		end = line[2]
		cStart = chrom + ":" + start
		cEnd = chrom + ":" + end
		if cStart in starts and cEnd in ends:
			line.append(starts[cStart])
			outJunc.append(line)

with open( sys.argv[3], 'w' ) as outFile:
    for x in outJunc:
        x = map(str, x)
        outLine = "\t".join(x)
        outFile.write( "%s\n" % ( outLine ) )
