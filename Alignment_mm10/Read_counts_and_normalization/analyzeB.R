
library(edgeR)
library(WGCNA)
library(biomaRt)
library(plotrix)

library(foreach)
library(doMC)
registerDoMC()
library(proxy)

getDoParWorkers()
options(cores=14)
getDoParWorkers()

setwd("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization")

source("functionDefinitions.R")
try(dir.create("resultsCoexpr"), silent = T)
try( dir.create("figuresCoexpr"), silent = F)
try(dir.create("resultsCoSplicEx"), silent = T)
try( dir.create("figuresCoSplicEx"), silent = F)

geneReadsRaw=as.matrix(read.table("data/RNA160225TP_gene_reads_not_normalized.txt", sep="\t",header=T, row.names=1))

sample_names<-read.table("data/samples.txt", sep=",")

# divide the data in different groups
H_group<-sample_names[,grepl("H",t(sample_names))]
L_group<-sample_names[,!grepl("H",t(sample_names))]
H_genes_counts<-geneReadsRaw[,colnames(geneReadsRaw) %in% t(H_group)]
L_genes_counts<-geneReadsRaw[,colnames(geneReadsRaw) %in% t(L_group)]
save(H_genes_counts,H_group,L_genes_counts,L_group, file="data/groups.RData")

groupSelection=c(rep("H_group",dim(H_genes_counts)[2]),rep("L_group",dim(L_genes_counts)[2]))
groupSelection =factor(groupSelection)

d=DGEList(counts= cbind(H_genes_counts, L_genes_counts), group= groupSelection)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")
sum(de.calls)
resultsDEtotal=cbind(de.tgw$table, de.calls)

# select genes with logCPM > 0 for further inclusion in network construction
resultsDEtotal=resultsDEtotal[resultsDEtotal[,"logCPM"]>0,]
write.csv(resultsDEtotal, file="resultsCoexpr/resultsDEtotal_083016.csv")
resultsDEtotal=read.csv("resultsCoexpr/resultsDEtotal_083016.csv")
rownames(resultsDEtotal)=resultsDEtotal[,1]

###################################################################################

# calculate edgeR normalization factors and normalize the data - use all data not just selected
UQnormFactors=calcNormFactors(geneReadsRaw, method=c("upperquartile"))

effectiveLibrarySizes= UQnormFactors*colSums(geneReadsRaw)
meanEffLibSize=mean(effectiveLibrarySizes)
countNormFactor= meanEffLibSize/effectiveLibrarySizes

normalizedGeneCountsUQ=0* geneReadsRaw

for (sample in 1:dim(normalizedGeneCountsUQ)[2]){  
  normalizedGeneCountsUQ[,sample]= geneReadsRaw[, sample]* countNormFactor[sample]	
}

geneReadsNorm=normalizedGeneCountsUQ[rownames(resultsDEtotal), ]
##################################################
# # find genes with high number of reads

connCounts=softConnectivity(t(geneReadsNorm), power=6)

quantileConn=quantile(connCounts, seq(0, 1, 0.1))  
geneReadsHighConn=geneReadsNorm[connCounts>quantileConn[6],]
geneNamesHighConn=rownames(geneReadsHighConn)
###################################################################################
exonReadsRaw=read.table("data/RNA160225TP_exon_reads_not_normalized.txt", sep="\t",header=T, row.names=1)

# UQnormFactors=calcNormFactors(exonReadsRaw, method=c("upperquartile"))
# 
# effectiveLibrarySizes= UQnormFactors*colSums(exonReadsRaw)
# meanEffLibSize=mean(effectiveLibrarySizes)
# countNormFactor= meanEffLibSize/effectiveLibrarySizes


exonReadsNorm=0* exonReadsRaw
for (sample in 1:dim(exonReadsRaw)[2]){
  exonReadsNorm[,sample]= exonReadsRaw[, sample]* countNormFactor[sample]	
}
exonReadsNorm =round(exonReadsNorm)

# select exons from genes with at least 1 CPM
splitIDs=mapply(strsplit, rownames(exonReadsNorm), MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 4))

exonCountsHighCounts=exonReadsNorm[which(exonGeneName %in% rownames(resultsDEtotal)),]
exonGeneNamesHighCounts=exonGeneName[exonGeneName %in% rownames(resultsDEtotal)]

canberraListExons=foreach (geneName = rownames(resultsDEtotal), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
  colnames(exonDistMatrix)=colnames(exonCountsHighCounts)
  rownames(exonDistMatrix)=colnames(exonCountsHighCounts)
  exonDistMatrix
  
}
names(canberraListExons)=rownames(resultsDEtotal)

########################################################################################################
########################################################################################################

save(canberraListExons, file="data/canberraListExons_083016.RData")
load("data/canberraListExons_083016.RData")

nGenes=length(canberraListExons)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExons[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExons))
colnames(distData)=names(canberraListExons)

for(gene in names(canberraListExons)) {
  distData[,gene]=as.vector(as.dist(canberraListExons[[gene]]))
} 

adjCoSplicEx_large=adjacency(distData,power=6)


save(adjCoSplicEx_large, file="data/adjCoSplicEx_large_083016.RData")

load("data/adjCoSplicEx_large_083016.RData")
#just in case ...
adjCoSplicEx_large[is.na(adjCoSplicEx_large)]=0
diag(adjCoSplicEx_large)=1
colnames(adjCoSplicEx_large)=rownames(adjCoSplicEx_large)


connCoSplicEx=rowSums(adjCoSplicEx_large)
quantileConnExons=quantile(connCoSplicEx, probs = seq(0, 1, 0.1))  

geneNamesHighCoSplicExConn=names(canberraListExons)[connCoSplicEx>quantileConnExons[6]]
################################################################################################3
selectedGeneCounts=geneReadsHighConn

canberraListSelected=canberraListExons[geneNamesHighCoSplicExConn]
adjCoSplicEx=adjCoSplicEx_large[geneNamesHighCoSplicExConn,geneNamesHighCoSplicExConn]

exonGeneNameSelected=geneNamesHighCoSplicExConn
selectedExonCounts=exonReadsNorm[which(exonGeneName %in% geneNamesHighCoSplicExConn),]

########################################################################################################
canberraListExonsHigh=foreach (geneName = rownames(adjCoSplicEx), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
  
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),sample_info[,"Line"]=="H"]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
  #   colnames(exonDistMatrix)=exonColnames
  #   rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExonsHigh)=rownames(adjCoSplicEx)
nGenes=length(canberraListExonsHigh)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExonsHigh[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExonsHigh))
colnames(distData)=names(canberraListExonsHigh)

for(gene in names(canberraListExonsHigh)) {
  distData[,gene]=as.vector(as.dist(canberraListExonsHigh[[gene]]))
} 

adjCoSplicEx_High=adjacency(distData,power=6)
adjCoSplicEx_High[is.na(adjCoSplicEx_High)]=0

canberraListExonsLow=foreach (geneName = rownames(adjCoSplicEx), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
  #currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),sample_info[,"Line"]=="L"]  
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),unlist(L_group)]  
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
  #   colnames(exonDistMatrix)=exonColnames
  #   rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExonsLow)=rownames(adjCoSplicEx)
nGenes=length(canberraListExonsLow)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExonsLow[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExonsLow))
colnames(distData)=names(canberraListExonsLow)

for(gene in names(canberraListExonsLow)) {
  distData[,gene]=as.vector(as.dist(canberraListExonsLow[[gene]]))
} 

adjCoSplicEx_Low=adjacency(distData,power=6)
adjCoSplicEx_Low[is.na(adjCoSplicEx_Low)]=0


########################################################################################################
########################################################################################################




save(selectedGeneCounts, canberraListSelected,adjCoSplicEx,adjCoSplicEx_Low,adjCoSplicEx_High,selectedExonCounts, exonGeneNameSelected,canberraListExonsLow,adjCoSplicEx_Low,file="data/selectedData.RData")
load("data/selectedData.RData")


