library(foreach)
library(doMC)
registerDoMC()
library(WGCNA)
library(multtest)

getDoParWorkers()
options(cores=4)
getDoParWorkers()
library(sgof)
library(edgeR)

setwd("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/GeneNetworkAnalysis")
source("functionDefinitions.R")
try(dir.create("resultsCoSplicEx"), silent = T)
try(dir.create("resultsCoSplicEx/moduleGeneList"), silent = T)


load("../data/selectedCountData.RData")
geneNames=rownames(adjCoSplicEx)
############################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjCoSplicEx, powerVector = powers, verbose = 5, moreNetworkConcepts=T)

plotNetConstruction(sft)
quartz.save("figures/netConstructionCoSplicEx.tif", type="tif", bg="white", dpi=300)
quartz.save("figures/netConstructionCoSplicEx.jpg", type="jpg", bg="white")

softPowerCoSplicEx=4
adjCoSplicEx=adjCoSplicEx^softPowerCoSplicEx
adjCoSplicEx[is.na(adjCoSplicEx)]=0
hierADJCoSplicEx = hclust(as.dist(1-adjCoSplicEx),method="average");

# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoSplicEx=cutreeHybrid(dendro = hierADJCoSplicEx, distM=1-adjCoSplicEx, cutHeight = 0.999999999995, maxCoreScatter=0.999995, minClusterSize = 20, deepSplit = 4, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, 
                             pamStage = T, maxPamDist=1,pamRespectsDendro = F, useMedoids = FALSE,  respectSmallClusters = T, verbose = 4, indent = 0)

colorsCoSplicEx = labels2colors(hybridCoSplicEx$labels)
names(colorsCoSplicEx)=geneNames
table(colorsCoSplicEx)
length(table(colorsCoSplicEx))
modulesCoSplicEx=names(table(colorsCoSplicEx))
sum(colorsCoSplicEx=="grey")
modulesCoSplicEx=names(table(colorsCoSplicEx))

save(adjCoSplicEx, colorsCoSplicEx,modulesCoSplicEx, file="data/adjModulesCoSplicEx.RData")
load("data/adjModulesCoSplicEx.RData")

########################################################################################################
# save the results below for use with enrinchR
CoSplicExConn=intramodularConnectivity(adjCoSplicEx, colorsCoSplicEx, scaleByMax=T)
totalScaledConnectivity=CoSplicExConn[,"kTotal"]/max(CoSplicExConn[,"kTotal"])
CoSplicExConn=cbind(CoSplicExConn, totalScaledConnectivity)

for (module in modulesCoSplicEx){
  print(module)
  currModuleInfo=cbind(rownames(CoSplicExConn)[colorsCoSplicEx==module],CoSplicExConn[colorsCoSplicEx==module,"kWithin"])
  write.csv(currModuleInfo, file=paste("resultsCoSplicEx/moduleGeneList/module_", module, ".csv", sep=""), row.names=F, col.names=F)  
}
#############################################################################
# GO annotations

load("data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsCoSplicEx, transcriptInfoMouse, "CoSplicEx")

##############################################################################
neuronsList=read.csv("data/CahoyNeurons.csv", header=TRUE)
neuronsSymbols= neuronsList[,"Gene.Name"]

astrosList=read.csv("data/CahoyAstros.csv", header=TRUE)
astrosSymbols= astrosList[,"Gene.Name"]

oligosList=read.csv("data/CahoyOligos.csv", header=TRUE)
oligosSymbols= oligosList[,"Gene.Name"]

moduleEnrichmentNeurons = moduleEnrichment (colorsCoSplicEx, neuronsSymbols)
moduleEnrichmentAstros = moduleEnrichment (colorsCoSplicEx, astrosSymbols)
moduleEnrichmentOligos = moduleEnrichment (colorsCoSplicEx, oligosSymbols)

cellTypeEnrichment=round(cbind(moduleEnrichmentNeurons,moduleEnrichmentAstros, moduleEnrichmentOligos),4)
colnames(cellTypeEnrichment)=c("Neurons", "Astros", "Oligos")
rownames(cellTypeEnrichment)=modulesCoSplicEx

write.csv(cellTypeEnrichment, file="resultsCoSplicEx/cellTypeEnrich.csv", append=T)

#####################################################################################################3


# evaluate splicing significance

HSCC_H=selectedExonCounts[,samplesHigh ]
HSCC_L=selectedExonCounts[,samplesLow]

############################################################################################
# evaluate genes with high splicing variability - this includes differential splicing AND differential variability in splicing
names(groupSelection)=colnames(selectedExonCounts)
phenotypeVector=factor(groupSelection)
names(phenotypeVector)=colnames(selectedExonCounts)
diffVarSplicingResults = diffVarSplicing(canberraListSelected, phenotypeVector, nPerm=200, nCores=27) 

save(diffVarSplicingResults, file="resultsCoSplicEx/diffVarSplicing.RData")
load("resultsCoSplicEx/diffVarSplicing.RData")

############################################################################################
# save summary differential splicing and differential splicing variability
colnames(diffVarSplicingResults)=c("A", "p value", "adjusted p")


summaryResults=cbind(geneNames, colorsCoSplicEx[geneNames], round(cbind(CoSplicExConn[geneNames,c("kWithin")], diffVarSplicingResults[geneNames,2:3]),3))

colnames(summaryResults)[3]="modularConnectivity"
colnames(summaryResults)[2]="module"

write.csv(summaryResults, file="resultsCoSplicEx/summaryResults.csv")



############################################################################################
# find modules enriched in differentially spliced/variable genes

pValues=summaryResults[ ,"p value"]
names(pValues)=geneNames
dSplicedGenes=names(pValues)[pValues<0.01]
write.csv(summaryResults[dSplicedGenes,], file="resultsCoSplicEx/affectedGenesDS.csv")

networkGenes=names(colorsCoSplicEx)
annotateMouseGeneGroupGO (dSplicedGenes,networkGenes, transcriptInfoMouse, type="CoSplicEx")


##################################################################################################



DSModuleEnrich=moduleEnrichment(colorsCoSplicEx, dSplicedGenes)
dsAffectedModules=names(DSModuleEnrich[DSModuleEnrich<(0.05/length(modulesCoSplicEx))])

fileConnSummary<-file("resultsCoSplicEx/SummaryResultsCoSplicEx.txt",  open="w")
writeLines(paste("Modules enriched in DS genes ", dsAffectedModules), fileConnSummary)
close(fileConnSummary)
# ##################################################################

#############################################################################################################333
# now do DE for each exon separately
exonGeneNames=rownames(selectedExonCounts)
exonGeneNamesSplit=strsplit(exonGeneNames, "_")
exonGeneNamesIdent=exonGeneNames
for (exon in 1:length(exonGeneNamesSplit)){
  exonGeneNamesIdent[exon]=exonGeneNamesSplit[[exon]][1]
  
}


d=DGEList(counts= cbind(HSCC_H, HSCC_L), group= groupSelection)
d <- estimateCommonDisp(d)

d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")
sum(de.calls)

pValuesRaw=de.tgw$table[,"PValue"]
names(pValuesRaw)=exonGeneNames

sum(pValuesRaw<0.01)

#####################################################################

sortIndexes=sort.int(pValuesRaw, decreasing = F, index.return=T)$ix
sortedGeneNames=exonGeneNames[sortIndexes]

adjustedResults<-SGoF(u=pValuesRaw)


summary(adjustedResults)

sortedPvals=pValuesRaw[sortedGeneNames]
sortedAdjustedPvals=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals)=names(sortedPvals)


####################################################################
# adjpOut=mt.rawp2adjp(pValuesRaw, proc=c( "BH"))
# fdrDE=adjpOut$adjp[order(adjpOut$index),2]
# names(fdrDE)=exonGeneNames
# sum(fdrDE<0.01)
# sum(fdrDE<0.05)
# fdrDE[fdrDE<0.05]

summaryResultsExon=cbind(colorsCoSplicEx[exonGeneNamesIdent],exonGeneNamesIdent, exonGeneNames, round(cbind(CoSplicExConn[exonGeneNamesIdent, "kWithin"], sortedPvals[exonGeneNames], sortedAdjustedPvals[exonGeneNames]), 3))
colnames(summaryResultsExon)[1]="CoSplicEx module"
colnames(summaryResultsExon)[4:6]=c("module conn", "p value", "adusted p")

summaryResultsExon=summaryResultsExon[summaryResultsExon[,4] > 0.7 & summaryResultsExon[,5] < 0.01,]

colnames(summaryResultsExon)[4:6]=c("connConsensus",  "P value", "adjusted P")
write.csv(summaryResultsExon, file="resultsCoSplicEx/summaryResultsExon.csv")

#############################################################################################################333
# now do DE for each exon separately