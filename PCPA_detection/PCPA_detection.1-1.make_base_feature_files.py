import sys
import re
import datetime

# Argument 1: Sample Name
# Argument 2: quartered intron/UTR coverage file from bedtools
# Argument 3: all feature coverage file from bedtools
# Argument 4: reference GFF3 file (refGene_full.gff3 - all genes, not longest isoform)
# Argument 5: STAR junction file (SJ.out.tab)
# Argument 6: GFold total mapped read normalized FPKM file (.read.count)
# Argument 7: Select feature lengths file for the genome of use (XYZ_select_feature_lengths)
# Argument 8: Mapped reads, per million for normalization

sample = sys.argv[1]
MMR = float(sys.argv[8])

# take the reference GFF3 file and convert it into a dict
print "Beginning analysis at %s" % (datetime.datetime.now())
geneRef = {} # make the dict for the gene reference
with open( sys.argv[4] ) as gff3Ref:
	next(gff3Ref)
	for i in gff3Ref:
		line = i.strip()
		line = line.split("\t")
		region = line[2]
		chrom = line[0]
		start = int(line[3])
		end = int(line[4])
		strand = line[6]
		length = end - start
		if region == "mRNA" or region == "transcript": # we're only interested in the genes themselves, not introns, exons, etc.
			parent = line[8].split("=")[2] # each isoform has a gene name and 
			geneID = line[8].split(";")[0].split("=")[1] # an isoform ID
			if parent not in geneRef:
				geneRef[parent] = {} # we need to create a dict for each gene to include feature info
				geneRef[parent]["geneLength"] = length
				geneRef[parent]["geneStart"] = start
				geneRef[parent]["geneEnd"] = end
				geneRef[parent]["geneStrand"] = strand
				geneRef[parent]["geneChrom"] = chrom
				geneRef[parent]["geneIDs"] = [geneID]
				geneRef[parent]["annotatedJuncReads"] = 0
				geneRef[parent]["novelJuncReads"] = 0
			elif parent in geneRef and length > geneRef[parent]["geneLength"]: # we only want to use info for the longest isoform
				geneRef[parent]["geneLength"] = length
				geneRef[parent]["geneStart"] = start
				geneRef[parent]["geneEnd"] = end
				geneRef[parent]["geneIDs"].append(geneID) # we keep each isoform ID for the genes
			else:
				geneRef[parent]["geneIDs"].append(geneID)

#### Begin quartered features section ####
print "Reading quartered features BED file at %s" % (datetime.datetime.now())
quarteredIntrons = {}
quarteredUTRs = {}
with open( sys.argv[2] ) as quarteredFeatures:
	for i in quarteredFeatures:
		line = i.strip()
		line = line.split("\t")
		region = line[2]
		geneName = line[8].split(":")[1]
		quarteredIntrons[geneName] = quarteredIntrons.get(geneName, {}) # if there isn't a dict already for this gene, create it
		quarteredUTRs[geneName] = quarteredUTRs.get(geneName, {}) # if there isn't a dict already for this gene, create it
		readCount = int(line[9])
		if region == "intron":
			intronSection = float(line[8].split(":")[2].split(";")[0]) # multiply the quartered feature by 100 to convert it to an int from a float
			quarteredIntrons[geneName][intronSection] = readCount
		elif region == "three_prime_UTR":
			utrSection = float(line[8].split(":")[2].split(";")[0])
			quarteredUTRs[geneName][utrSection] = readCount

for gene in quarteredUTRs: # if the UTR is split over 2 exons, combine them into one UTR
	if 2.0 in quarteredUTRs[gene]:
		quarteredUTRs[gene][1.0] = quarteredUTRs[gene][1.0] + quarteredUTRs[gene][1.25] # combine congruous sections
		quarteredUTRs[gene][1.25] = quarteredUTRs[gene][1.5] + quarteredUTRs[gene][1.75]
		quarteredUTRs[gene][1.5] = quarteredUTRs[gene][2.0] + quarteredUTRs[gene][2.25]
		quarteredUTRs[gene][1.75] = quarteredUTRs[gene][2.5] + quarteredUTRs[gene][2.75]
		del quarteredUTRs[gene][2.0] # remove the old key:value pairs
		del quarteredUTRs[gene][2.25]
		del quarteredUTRs[gene][2.5]
		del quarteredUTRs[gene][2.75]
			
for gene in geneRef:
	for isoform in geneRef[gene]["geneIDs"]:
		if isoform in quarteredIntrons:
			geneRef[gene]["intron1Quarter1Reads"] = quarteredIntrons[isoform].get("1.0","#N/A")
			geneRef[gene]["intron1Quarter2Reads"] = quarteredIntrons[isoform].get("1.25","#N/A")
			geneRef[gene]["intron1Quarter3Reads"] = quarteredIntrons[isoform].get("1.5","#N/A")
			geneRef[gene]["intron1Quarter4Reads"] = quarteredIntrons[isoform].get("1.75","#N/A")
			geneRef[gene]["intron2Quarter1Reads"] = quarteredIntrons[isoform].get("2.0","#N/A")
			geneRef[gene]["intron2Quarter2Reads"] = quarteredIntrons[isoform].get("2.25","#N/A")
			geneRef[gene]["intron2Quarter3Reads"] = quarteredIntrons[isoform].get("2.5","#N/A")
			geneRef[gene]["intron2Quarter4Reads"] = quarteredIntrons[isoform].get("2.75","#N/A")
		if isoform in quarteredUTRs:
			geneRef[gene]["3UTRQuarter1Reads"] = quarteredUTRs[isoform].get("1.0","#N/A")
			geneRef[gene]["3UTRQuarter2Reads"] = quarteredUTRs[isoform].get("1.25","#N/A")
			geneRef[gene]["3UTRQuarter3Reads"] = quarteredUTRs[isoform].get("1.5","#N/A")
			geneRef[gene]["3UTRQuarter4Reads"] = quarteredUTRs[isoform].get("1.75","#N/A")
#### End quartered features section ####


#### Begin all features section ####
# create the dictionaries used for the read values
print "Reading all features BED file at %s" % (datetime.datetime.now())
genesIntrons = {}
genesExons = {}
genesCDS = {}
genes3UTRs = {}
genes5UTRs = {}
genesTranscript = {}
with open( sys.argv[3] ) as allFeatures:
	for i in allFeatures:
		line = i.strip()
		line = line.split("\t")
		region = line[2] # intron, exon, CDS, mRNA, transcript, UTRs
		geneName = (line[8].split(";")[0].split(":")[1] if line[8].__contains__(":") else line[8].split(";")[0].split("=")[1])  # changes the split to extract the gene ID based on the identifier field
		readCount = int(line[9])
		if region == "intron":
			genesIntrons[geneName] = genesIntrons.get(geneName, {}) # if there isn't a dict already for this gene, create it
			intronNum = int(line[8].split(";")[0].split(":")[2]) # extract intron number
			genesIntrons[geneName][intronNum] = readCount # nested dict: gene -> intron number -> read count
		elif region == "exon":
			genesExons[geneName] = genesExons.get(geneName, {}) # if there isn't a dict already for this gene, create it
			exonNum = int(line[8].split(";")[0].split(":")[2]) # extract exon number
			genesExons[geneName][exonNum] = readCount # nested dict: gene -> exon number -> read count
		elif region == "CDS":
			genesCDS[geneName] = genesCDS.get(geneName, {}) # if there isn't a dict already for this gene, create it
			cdsNum = int(line[8].split(";")[0].split(":")[2]) # extract CDS number
			genesCDS[geneName][cdsNum] = readCount # nested dict: gene -> CDS number -> read count
		elif region == "three_prime_UTR":
			genes3UTRs[geneName]  = genes3UTRs.get(geneName, 0) # if there isn't a value already for this gene, create it
			genes3UTRs[geneName] += readCount # add the new value to the old
		elif region == "five_prime_UTR":
			genes5UTRs[geneName] = genes5UTRs.get(geneName, 0) # if there isn't a value already for this gene, create it
			genes5UTRs[geneName] += readCount # add the new value to the old
		elif region == "mRNA" or region == "transcript":
			genesTranscript[geneName] = genesTranscript.get(geneName, 0) # if there isn't a value already for this gene, create it
			if readCount > genesTranscript[geneName]: # check for read count, if it's more than the current replace
				genesTranscript[geneName] = readCount
				
# now we check the genes against the feature reads
# this will append the relevant info to the geneRef dict
print "Appending feature information to each gene at %s" % (datetime.datetime.now())
for gene in geneRef:
	for isoform in geneRef[gene]["geneIDs"]: # we need to check each isoform against the features
		if isoform in genesTranscript:
			# these check if there is a feature (intron here) associated with the relevant number in the gene
			# for example, it checks the isoform ID in the intron dict and then looks for intron 1
			# if intron 1 isn't found, it adds #N/A as the alternative
			geneRef[gene]["firstIntronReads"] = genesIntrons.get(isoform, {}).get(1,"#N/A")
			geneRef[gene]["secondIntronReads"] = genesIntrons.get(isoform, {}).get(2,"#N/A")
			geneRef[gene]["thirdIntronReads"] = genesIntrons.get(isoform, {}).get(3,"#N/A")
			try:
				geneRef[gene]["totalIntronReads"] =  sum(genesIntrons[isoform][intron] for intron in genesIntrons[isoform])
			except KeyError:
				geneRef[gene]["totalIntronReads"] =  "#N/A"
			# this checks for the feature associated with the isoform ID
			# if they are found, it sums them up for the total
			# if not (KeyError) it just adds #N/A
			# same thing for exon
			geneRef[gene]["firstExonReads"] = genesExons.get(isoform, {}).get(1,"#N/A")
			geneRef[gene]["secondExonReads"] = genesExons.get(isoform, {}).get(2,"#N/A")
			geneRef[gene]["thirdExonReads"] = genesExons.get(isoform, {}).get(3,"#N/A")
			try:
				geneRef[gene]["totalExonReads"] =  sum(genesExons[isoform][exon] for exon in genesExons[isoform])
			except KeyError:
				geneRef[gene]["totalExonReads"] =  "#N/A"
			# this checks for the last feature based on key number in the feature dict
			# if there is a max (i.e. last) it counts that as the terminal feature
			try:
				geneRef[gene]["lastExonReads"] = max([exon for exon in genesExons[isoform].values()])
			except KeyError:
				geneRef[gene]["lastExonReads"] = "#N/A"
			# same as above for CDS
			geneRef[gene]["firstCDSReads"] = genesCDS.get(isoform, {}).get(1,"#N/A")
			geneRef[gene]["secondCDSReads"] = genesCDS.get(isoform, {}).get(2,"#N/A")
			geneRef[gene]["thirdCDSReads"] = genesCDS.get(isoform, {}).get(3,"#N/A")
			try:
				geneRef[gene]["totalCDSReads"] =  sum(genesCDS[isoform][cds] for cds in genesCDS[isoform])
			except KeyError:
				geneRef[gene]["totalCDSReads"] =  "#N/A"
			try:
				geneRef[gene]["lastCDSReads"] = max([cds for cds in genesCDS[isoform].values()])
			except KeyError:
				geneRef[gene]["lastCDSReads"] = "#N/A"
			# we already have summed the 3' and 5' UTRs as well as all gene reads
			geneRef[gene]["3UTRReads"] = genes3UTRs.get(isoform,"#N/A")
			geneRef[gene]["5UTRReads"] = genes5UTRs.get(isoform,"#N/A")
			geneRef[gene]["totalReads"] = genesTranscript.get(isoform,"#N/A")

#### End all features section ####

#### Begin junction section ####
# take the STAR junction file and convert it into a dict
print "Reading STAR junction file at %s" % (datetime.datetime.now())
spliceJunc = {}
with open( sys.argv[5] ) as starJunctions: # I converted to using the STAR junction files
	juncCount = 1
	for i in starJunctions:
		line = i.strip()
		line = line.split("\t")
		chrom = line[0]
		start = int(line[1])
		end = int(line[2])
		strand = line[3]
		notation = line[5]
		juncReads = int(line[6]) + int(line[7])
		spliceJunc[juncCount] = {}
		spliceJunc[juncCount]["chrom"] = chrom
		spliceJunc[juncCount]["start"] = start
		spliceJunc[juncCount]["end"] = end
		spliceJunc[juncCount]["strand"] = "+" if strand is "1" else "-"
		spliceJunc[juncCount]["notation"] = "unannotated" if notation is "0" else "annotated" # convert STAR notation into string
		spliceJunc[juncCount]["reads"] = juncReads
		juncCount += 1

for junc in spliceJunc: # this is a long double loop
	for gene in geneRef: # I can't find a way around it
		if spliceJunc[junc]["chrom"] == geneRef[gene]["geneChrom"]: # it checks each spliced junction location against a gece location
			if spliceJunc[junc]["start"] >= geneRef[gene]["geneStart"] and spliceJunc[junc]["end"] <= geneRef[gene]["geneEnd"]: # have to check the junction location against gene start/ends
				if spliceJunc[junc]["notation"] == "annotated":
					geneRef[gene]["annotatedJuncReads"] += spliceJunc[junc]["reads"]
				else:
					geneRef[gene]["novelJuncReads"] += spliceJunc[junc]["reads"]
					
""" Unused as of now, I converted to the STAR Junction files
exonExonJuncs = {}
with open( sys.argv[4] ) as exonJunctions:
	for i in exonJunctions:
		line = i.strip()
		line = line.split("\t")
		if line[21] is not "0":
			juncID = line[3]
			juncReads = int(line[4])
			geneID = line[20].split(":")[1]
			exonNum = line[20].split(":")[2].split(";")[0]
			if juncID not in exonExonJuncs

geneJuncs = {}
with open( sys.argv[5] ) as geneJunctions:
	for i in geneJunctions:
		line = i.strip()
		line = line.split("\t")
		juncID = line[3]
		juncReads = int(line[4])
		geneID = line[20].split(":")[1]
		exonNum = line[20].split(":")[2].split(";")[0]
"""
#### End junction section ####


#### Begin FPKM section ####
print "Reading GFold FPKM count file at %s" % (datetime.datetime.now())
with open( sys.argv[6] ) as GFold:
	for i in GFold:
		line = i.strip()
		line = line.split("\t")
		gene = line[0]
		fpkm = line[4]
		if gene in geneRef:
			geneRef[gene]["FPKM"] = fpkm
#### End FPKM section ####

#### Add in feature lengths ####
print "Adding gene/feature lenght information at %s" % (datetime.datetime.now())
geneFeatures = {}
with open( sys.argv[7] ) as lengthData:
	next(lengthData)
	for i in lengthData:
		line = i.strip()
		line = line.split("\t")
		geneFeatures[line[0]] = {} # make a dict for each gene ID
		geneFeatures["mRNALength"] = line[2]
		geneFeatures["exon1Length"] = line[3]
		geneFeatures["exon2Length"] = line[4]
		geneFeatures["exon3Length"] = line[5]
		geneFeatures["intron1Length"] = line[6]
		geneFeatures["intron2Length"] = line[7]

for gene in geneRef:
	for isoform in geneRef[gene]["geneIDs"]:
		if isoform in geneFeatures:
			geneRef[gene]["mRNALength"] = geneFeatures[isoform].get("mRNALength","#N/A")
			geneRef[gene]["exon1Length"] = geneFeatures[isoform].get("exon1Length","#N/A")
			geneRef[gene]["exon2Length"] = geneFeatures[isoform].get("exon2Length","#N/A")
			geneRef[gene]["exon3Length"] = geneFeatures[isoform].get("exon3Length","#N/A")
			geneRef[gene]["intron1Length"] = geneFeatures[isoform].get("intron1Length","#N/A")
			geneRef[gene]["intron2Length"] = geneFeatures[isoform].get("intron2Length","#N/A")
#### End feature lengths ####

#### Put everything together ####
print "Printing output at %s" % (datetime.datetime.now())
with open( sample + ".raw_feature_file.txt", 'w' ) as outFile:
	outFile.write( "gene\tchrom\tgeneLength\tmRNALength\texon1Length\texon2Length\texon3Length\tintron1Length\tintron2Length\ttotalGeneReads\
\tFPKM\tannotatedJunctions\tnovelJunctions\tfirstCDS\tsecondCDS\tthirdCDS\tlastCDS\ttotalCDS\tfirstExon\tsecondExon\tthirdExon\tlastExon\ttotalExon\
\tfirstIntron\tsecondIntron\tthirdIntron\ttotalIntron\t5UTR\t3UTR\tintron1Quarter1\tintron1Quarter2\tintron1Quarter3\tintron1Quarter4\tintron2Quarter1\
\tintron2Quarter2\tintron2Quarter3\tintron2Quarter4\t3UTRQuarter1\t3UTRQuarter2\t3UTRQuarter3\t3UTRQuarter4\n" )
	for gene in geneRef:
		outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, geneRef[gene]["geneChrom"], geneRef[gene]["geneLength"], geneRef[gene]["mRNALength"], geneRef[gene]["exon1Length"], geneRef[gene]["exon2Length"], geneRef[gene]["exon3Length"], geneRef[gene]["intron1Length"], geneRef[gene]["intron2Length"], geneRef[gene]["totalReads"], geneRef[gene]["FPKM"], geneRef[gene]["annotatedJuncReads"], geneRef[gene]["novelJuncReads"], geneRef[gene]["firstCDSReads"], geneRef[gene]["secondCDSReads"], geneRef[gene]["thirdCDSReads"], geneRef[gene]["lastCDSReads"], geneRef[gene]["totalCDSReads"], geneRef[gene]["firstExonReads"], geneRef[gene]["secondExonReads"], geneRef[gene]["thirdExonReads"], geneRef[gene]["lastExonReads"],geneRef[gene]["totalExonReads"],  geneRef[gene]["firstIntronReads"], geneRef[gene]["secondIntronReads"], geneRef[gene]["thirdIntronReads"], geneRef[gene]["totalIntronReads"], geneRef[gene]["5UTRReads"], geneRef[gene]["3UTRReads"], geneRef[gene].get("intron1Quarter1Reads","#N/A"), geneRef[gene].get("intron1Quarter2Reads","#N/A"), geneRef[gene].get("intron1Quarter3Reads","#N/A"), geneRef[gene].get("intron1Quarter4Reads","#N/A"), geneRef[gene].get("intron2Quarter1Reads","#N/A"), geneRef[gene].get("intron2Quarter2Reads","#N/A"), geneRef[gene].get("intron2Quarter3Reads","#N/A"), geneRef[gene].get("intron2Quarter4Reads","#N/A"), geneRef[gene].get("3UTRQuarter1Reads","#N/A"), geneRef[gene].get("3UTRQuarter2Reads","#N/A"), geneRef[gene].get("3UTRQuarter3Reads","#N/A"), geneRef[gene].get("3UTRQuarter4Reads","#N/A")))

# normalized to mapped reads per million
for gene in geneRef:
	for key in geneRef[gene]:
		if key.__contains__("Quarter"):
			try:
				geneRef[gene][key] = ("#N/A" if geneRef[gene][key] == "#N/A" else float(geneRef[gene][key])/MMR)
			except KeyError:
				geneRef[gene][key] = "#N/A"
		elif key.__contains__("Reads"):
			geneRef[gene][key] = ("#N/A" if geneRef[gene][key] == "#N/A" else float(geneRef[gene][key])/MMR)

with open( sample + ".normalized_feature_file.txt", 'w' ) as outFile:
	outFile.write( "gene\tchrom\tgeneLength\tmRNALength\texon1Length\texon2Length\texon3Length\tintron1Length\tintron2Length\ttotalGeneReads\
\tFPKM\tannotatedJunctions\tnovelJunctions\tfirstCDS\tsecondCDS\tthirdCDS\tlastCDS\ttotalCDS\tfirstExon\tsecondExon\tthirdExon\tlastExon\ttotalExon\
\tfirstIntron\tsecondIntron\tthirdIntron\ttotalIntron\t5UTR\t3UTR\tintron1Quarter1\tintron1Quarter2\tintron1Quarter3\tintron1Quarter4\tintron2Quarter1\
\tintron2Quarter2\tintron2Quarter3\tintron2Quarter4\t3UTRQuarter1\t3UTRQuarter2\t3UTRQuarter3\t3UTRQuarter4\n" )
	for gene in geneRef:
		outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, geneRef[gene]["geneChrom"], geneRef[gene]["geneLength"], geneRef[gene]["mRNALength"], geneRef[gene]["exon1Length"], geneRef[gene]["exon2Length"], geneRef[gene]["exon3Length"], geneRef[gene]["intron1Length"], geneRef[gene]["intron2Length"], geneRef[gene]["totalReads"], geneRef[gene]["FPKM"], geneRef[gene]["annotatedJuncReads"], geneRef[gene]["novelJuncReads"], geneRef[gene]["firstCDSReads"], geneRef[gene]["secondCDSReads"], geneRef[gene]["thirdCDSReads"], geneRef[gene]["lastCDSReads"], geneRef[gene]["totalCDSReads"], geneRef[gene]["firstExonReads"], geneRef[gene]["secondExonReads"], geneRef[gene]["thirdExonReads"], geneRef[gene]["lastExonReads"],geneRef[gene]["totalExonReads"],  geneRef[gene]["firstIntronReads"], geneRef[gene]["secondIntronReads"], geneRef[gene]["thirdIntronReads"], geneRef[gene]["totalIntronReads"], geneRef[gene]["5UTRReads"], geneRef[gene]["3UTRReads"], geneRef[gene].get("intron1Quarter1Reads","#N/A"), geneRef[gene].get("intron1Quarter2Reads","#N/A"), geneRef[gene].get("intron1Quarter3Reads","#N/A"), geneRef[gene].get("intron1Quarter4Reads","#N/A"), geneRef[gene].get("intron2Quarter1Reads","#N/A"), geneRef[gene].get("intron2Quarter2Reads","#N/A"), geneRef[gene].get("intron2Quarter3Reads","#N/A"), geneRef[gene].get("intron2Quarter4Reads","#N/A"), geneRef[gene].get("3UTRQuarter1Reads","#N/A"), geneRef[gene].get("3UTRQuarter2Reads","#N/A"), geneRef[gene].get("3UTRQuarter3Reads","#N/A"), geneRef[gene].get("3UTRQuarter4Reads","#N/A")))
#### Done ####
