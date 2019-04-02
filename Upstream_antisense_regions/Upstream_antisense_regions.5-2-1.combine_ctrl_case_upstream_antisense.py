import re
import sys

PCPAID = []
with open( sys.argv[1] ) as PCPAFile:
    for i in PCPAFile:
        line = i.strip()
        PCPAID.append(line)

upID = []
with open( sys.argv[2] ) as upFile:
    for i in upFile:
        line = i.strip()
        upID.append(line)

downID = []
with open( sys.argv[3] ) as downFile:
    for i in downFile:
        line = i.strip()
        downID.append(line)


with open( sys.argv[4] ) as featureFile:
    for i in featureFile:
        line = i.strip()
        line = line.split("\t")
        geneID = line[0]
        if geneID == "gene_ID":
            UAchange = "upstream_antisense_change"
            status = "gene_status"
            line.append(UAchange)
            line.append(status)
            print('\t'.join(map(str,line)))
        else:
            ctrlUA = line[12]
            caseUA = line[21]
            ctrlFPKM = float(line[5])
            caseFPKM = float(line[14])
            if ctrlUA == "#N/A":
                UAchange = "#N/A"
            elif ctrlUA == "0":
                UAchange = caseUA
            else:
                UAchange = float(caseUA) / float(ctrlUA)
            if geneID in PCPAID:
                status = "PCPAed"
            elif geneID in upID:
                status = "upregulated"
            elif geneID in downID:
                status = "downregulated"
            elif (ctrlFPKM  >=1.0) or (caseFPKM >=1.0):
                status = "expressed_no_change"
            else:
                status = "not_expressed"
            line.append(UAchange)
            line.append(status)
            print('\t'.join(map(str,line)))

