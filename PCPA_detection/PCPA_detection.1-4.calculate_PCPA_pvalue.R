## Calculating PCPA genes, from feature files
## made with help from Larry Singh and Chao Di

args <- commandArgs(TRUE)
if(length(args)==0){
	print("No files provided as arguments!")
	quit()
}

require(exactci)
require(plyr)
require(parallel)

# read in PCPA feature files
featureFile_ctrl <- read.table(args[1],head=T)
featureFile_u1amo <- read.table(args[2],head=T)

prefix <- featureFile_ctrl[,1:5]
ctrl_feature <- featureFile_ctrl[,6:dim(featureFile_ctrl)[2]]
u1amo_feature <- featureFile_u1amo[,6:dim(featureFile_u1amo)[2]]
colnames(ctrl_feature) <- paste('ctrl',colnames(ctrl_feature),sep="_")
colnames(u1amo_feature) <- paste('u1amo',colnames(u1amo_feature),sep="_")

## combine two files
ctrlVSamo <- cbind(prefix,ctrl_feature,u1amo_feature)

## step0: add 1 to reads count in case there is 0 value, add pseudocount 0.001 to normalized values
ctrlVSamo[,c("ctrl_quarter1","ctrl_quarter4","ctrl_lastExon","u1amo_quarter1","u1amo_quarter4","u1amo_lastExon")] = ctrlVSamo[,c("ctrl_quarter1","ctrl_quarter4","ctrl_lastExon","u1amo_quarter1","u1amo_quarter4","u1amo_lastExon")] + 1
ctrlVSamo[,c("ctrl_norm_quarter1","ctrl_norm_quarter4","ctrl_norm_lastExon","u1amo_norm_quarter1","u1amo_norm_quarter4","u1amo_norm_lastExon")] = ctrlVSamo[,c("ctrl_norm_quarter1","ctrl_norm_quarter4","ctrl_norm_lastExon","u1amo_norm_quarter1","u1amo_norm_quarter4","u1amo_norm_lastExon")] + 0.01

## step1: add columns of divisions 
ctrlVSamo$ctrl_intron_Q1_Q4 = ctrlVSamo$ctrl_quarter1 / ctrlVSamo$ctrl_quarter4
ctrlVSamo$u1amo_intron_Q1_Q4 = ctrlVSamo$u1amo_quarter1 / ctrlVSamo$u1amo_quarter4
# intron Q1/Q4 ratio of u1amo to ctrl
ctrlVSamo$u1amo_Q1Q4_VS_ctrl_Q1Q4 = ctrlVSamo$u1amo_intron_Q1_Q4 / ctrlVSamo$ctrl_intron_Q1_Q4
# norm_nextExon ratio of ctrl to u1amo
ctrlVSamo$ctrl_lastExon_VS_u1amo_lastExon = ctrlVSamo$ctrl_norm_lastExon / ctrlVSamo$u1amo_norm_lastExon


## step2: rpkm >= 1 in either ctrl or u1amo
# new:rpkm>=1 only in control
ctrlVSamo_filter1 <- ctrlVSamo[which(ctrlVSamo$ctrl_rpkm >= 1),] 

print("rpkm>=1 number:")
length(unique(ctrlVSamo_filter1[,1]))

## step3: more than 10 reads on the 1st quarter of intron 1 in u1amo, and 1st quarter of intron > 4th quarter of intron in U1 AMO
ctrlVSamo_filter2 <- ctrlVSamo_filter1[which(ctrlVSamo_filter1$u1amo_quarter1 >= 10  & ctrlVSamo_filter1$u1amo_quarter1 > ctrlVSamo_filter1$u1amo_quarter4),]
length(unique(ctrlVSamo_filter2[,1]))

## step4: statistic tests for the ratio/proportion of intron_1stQuarter to intron_4thQuarter in ctrl and u1amo
# Function for computing the statistics for each row
# cq1, cq4, cne, uq1, uq4, une are column names whose reads are going to test, i.e., "ctrl_quarter1","ctrl_quarter4","ctrl_nextExon","u1amo_quarter1","u1amo_quarter4","u1amo_nextExon" 
ctrl_total_reads <- as.integer(args[3])
u1amo_total_reads <- as.integer(args[4])

getStats <- function(i, x, cne, une, cq1, cq4, uq1, uq4) {
         poissonP <- poisson.exact(c(x[i,cne], x[i,une]),c(ctrl_total_reads,u1amo_total_reads),alternative="greater")$p.value
        # use fisher-exact test for 1stQ increase relative to 4thQ
        # Build a 2x2 contingency table
        countTable <- matrix(c(x[i,uq1],x[i,cq1],x[i,uq4],x[i,cq4]),2,2)
	fisherP <- fisher.test(countTable,alternative="greater")$p.value

        data.frame(poissonP=poissonP,fisherP=fisherP)
}
## Function for Multiple testing
multiTest <- function(x) {
        poissonQ <- p.adjust(x[,"poissonP"], method="BH")
        fisherQ <- p.adjust(x[,"fisherP"], method="BH")
        data.frame(x, poissonQ=poissonQ, fisherQ=fisherQ)
}



# get pvalues and qvalues for one file
ctrlVSamo.p <- do.call(rbind, mclapply(1:nrow(ctrlVSamo_filter2), getStats, x=ctrlVSamo_filter2, "ctrl_lastExon","u1amo_lastExon","ctrl_quarter1","ctrl_quarter4","u1amo_quarter1","u1amo_quarter4",  mc.cores=20)) # number of cores to use
# qvalue, Multiple testing
ctrlVSamo.q <- cbind(ctrlVSamo_filter2,multiTest(ctrlVSamo.p))


# step5: p-value or q-value <= cutoff
pcpa_gene <- ctrlVSamo.q[which(ctrlVSamo.q$poissonQ <= 0.05 & ctrlVSamo.q$fisherQ <=0.05),]
length(unique(pcpa_gene[,1]))


write.table(pcpa_gene, file=sub("\\.feature","\\.pcpa.txt",args[2]), sep="\t", quote=F, row.names=FALSE)

