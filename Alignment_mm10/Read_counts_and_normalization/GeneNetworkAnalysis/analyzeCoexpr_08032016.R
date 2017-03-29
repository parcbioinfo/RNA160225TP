library(foreach)
library(doMC)
registerDoMC()
library(multtest)
library(WGCNA)
library("org.Mm.eg.db")
library(biomaRt)
library(GOstats)
library("org.Mm.eg.db")
library("edgeR")
library(vegan)
library(sgof)

library(lawstat)

#source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
#biocLite("GOstats")

getDoParWorkers()
options(cores=4)
getDoParWorkers()

enableWGCNAThreads(nThreads = 4)

setwd("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/GeneNetworkAnalysis")
load("../data/ReadCountsNormalization.RData")
load("../data/selectedCountData.RData")
try(dir.create("resultsCoexpr"), silent = T)
try( dir.create("figuresCoexpr"), silent = F)


geneNames=rownames(normalizedGeneCountsUQ_genes_extra_selected)

####################################################################################33

H_genes_counts<-normalizedGeneCountsUQ_genes_extra_selected[,colnames(normalizedGeneCountsUQ_genes_extra_selected) %in% t(H_group)]
L_genes_counts<-normalizedGeneCountsUQ_genes_extra_selected[,colnames(normalizedGeneCountsUQ_genes_extra_selected) %in% t(L_group)]

########################################################################################################################################
adjConsensus=adjacency(t(normalizedGeneCountsUQ_genes_extra_selected), power=1)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjConsensus, powerVector = powers, verbose = 5, moreNetworkConcepts=T)

plotNetConstruction(sft)
quartz.save("figuresCoexpr/netConstructionCoexpr.tif", type="tif", bg="white", dpi=300)
quartz.save("figuresCoexpr/netConstructionCoexpr.jpg", type="jpg", bg="white")

softPowerCoexpr=7
adjCoexpr=adjConsensus^softPowerCoexpr
adjCoexpr[is.na(adjCoexpr)]=0
hierADJCoexpr = hclust(as.dist(1-adjCoexpr),method="average");

# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoexpr=cutreeHybrid(dendro = hierADJCoexpr, distM=1-adjCoexpr, cutHeight = 0.9995, minClusterSize = 100, deepSplit = 4, maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, pamStage = TRUE, pamRespectsDendro = F, useMedoids = FALSE,  respectSmallClusters = TRUE, verbose = 2, indent = 0)

colorsCoexpr = labels2colors(hybridCoexpr$labels)
names(colorsCoexpr)=geneNames
table(colorsCoexpr)
length(table(colorsCoexpr))
modulesCoexpr=names(table(colorsCoexpr))
sum(colorsCoexpr=="grey")

coexprConnConsensus=intramodularConnectivity(adjCoexpr,  colorsCoexpr, scaleByMax=T)
totalScaledConnectivity=coexprConnConsensus[,"kTotal"]/max(coexprConnConsensus[,"kTotal"])



fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Number modules ", length(table(colorsCoexpr)), sep=','), fileConnSummary)
writeLines(paste("Number grey genes  ",  sum(colorsCoexpr=="grey"), sep=','), fileConnSummary)

close(fileConnSummary)

names(colorsCoexpr)=geneNames
save(adjCoexpr, colorsCoexpr,  file="../data/adjCoexprModules.RData")
load("../data/adjCoexprModules.RData")


########################################################################################################
########################################################################################################
########################################################################################################
# save the results below for use with enrinchR

#############################################################################
# GO annotations

load("../data/transcriptInfoMouse.RData")
source("functionDefinitions.R")
annotateMouseModulesGO(colorsCoexpr, transcriptInfoMouse, type="Coexpr")

##############################################################################
########################################################################################################
# overlap between modules and neuronal genes list

neuronsList=read.csv("../data/CahoyNeurons.csv", header=TRUE)
neuronsSymbols= neuronsList[,"Gene.Name"]

astrosList=read.csv("../data/CahoyAstros.csv", header=TRUE)
astrosSymbols= astrosList[,"Gene.Name"]

oligosList=read.csv("../data/CahoyOligos.csv", header=TRUE)
oligosSymbols= oligosList[,"Gene.Name"]

moduleEnrichmentNeurons = moduleEnrichment (colorsCoexpr, neuronsSymbols)
moduleEnrichmentAstros = moduleEnrichment (colorsCoexpr, astrosSymbols)
moduleEnrichmentOligos = moduleEnrichment (colorsCoexpr, oligosSymbols)

cellTypeEnrichment=round(cbind(moduleEnrichmentNeurons,moduleEnrichmentAstros, moduleEnrichmentOligos),4)
colnames(cellTypeEnrichment)=c("Neurons", "Astros", "Oligos")
rownames(cellTypeEnrichment)=t(modulesCoexpr)


 fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
 writeLines(paste("Modules enriched in neuronal cell types\n "), fileConnSummary)
 close(fileConnSummary)

write.csv(cellTypeEnrichment, file="resultsCoexpr/cellTypeEnrich.csv", append=T)


########################################################################################################
# evaluate differential expression
groupSelection=c(rep("H_group",dim(H_genes_counts)[2]),rep("L_group",dim(L_genes_counts)[2]))
groupSelection =factor(groupSelection)


d=DGEList(counts= cbind(H_genes_counts, L_genes_counts), group= groupSelection)
d <- estimateCommonDisp(d)

d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")
sum(de.calls)

pValues=de.tgw$table[,"PValue"]
names(pValues)=geneNames

# adjpOut=mt.rawp2adjp(pValuesRaw, proc=c( "BH"))
# fdrDE=adjpOut$adjp[order(adjpOut$index),2]
# names(fdrDE)=geneNamesRaw

# results out of sgof come out sorted!! therefore raw p values also need to get sorted, and original gene names are put in sorted order


sortIndexes=sort.int(pValues, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNames[sortIndexes]

adjustedResults<-SGoF(u=pValues)

summary(adjustedResults)

sortedPvals=pValues[sortedGeneNames]
sortedAdjustedPvals=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals)=sortedGeneNames


meanH_genes_counts=rowMeans(H_genes_counts[geneNames,])
meanL_genes_counts=rowMeans(L_genes_counts[geneNames,])

summaryResultsDE=signif(cbind(coexprConnConsensus[,"kWithin"], sortedPvals[geneNames], sortedAdjustedPvals[geneNames]),3)

##############################################################################
# find differentially variable genes

pvalVar=rep(1, length(geneNames))
names(pvalVar)=geneNames

for (gene in geneNames){
  pvalVar[gene]=var.test(x=as.matrix(H_genes_counts[gene,]), y=as.matrix(L_genes_counts[gene,]))$p.value
}

sortIndexes=sort.int(pvalVar, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNames[sortIndexes]

adjustedResults<-SGoF(u=pvalVar)


summary(adjustedResults)

sortedPvals=pvalVar[sortedGeneNames]
sortedAdjustedPvals=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals)=names(sortedPvals)
# THIS IS NEVER USED
#sdH=apply(H_genes_counts, 1, sd)
#sdL=apply(L_genes_counts, 1, sd)

summaryResults=cbind(colorsCoexpr, round(summaryResultsDE, 3), round(cbind(sortedPvals[geneNames], sortedAdjustedPvals[geneNames]),3))
colnames(summaryResults)=c("colorsCoexpr","connConsensus", "p DE", "adjusted p DE", "p DV", "adjusted p DV")

write.csv(summaryResults, file="resultsCoexpr/SummaryResultsDEDV.csv")

##############################################################################
# evaluate module enrichment in differentially expressed/variable genes


#################
pValuesDE=summaryResults[ ,"p DE"]
names(pValuesDE)=geneNames
deGenes=names(pValuesDE)[pValuesDE<0.01]

pValuesDV=summaryResults[ ,"p DV"]
names(pValuesDV)=geneNames
dvGenes=names(pValuesDV)[pValuesDV<0.01]

affectedGenes=union(deGenes, dvGenes)
write.csv(summaryResults[affectedGenes,], file="resultsCoexpr/affectedGenesDEDV.csv")
affectedGenes=read.csv("resultsCoexpr/affectedGenesDEDV.csv")[,1]

networkGenes=names(colorsCoexpr)
annotateMouseGeneGroupGO(affectedGenes,networkGenes, transcriptInfoMouse, type="Coexpr")


#################

# ##############################################################################
##############################################################################

moduleEnrich=moduleEnrichment(colorsCoexpr, union(deGenes, dvGenes))
affectedModules=names(moduleEnrich[moduleEnrich<(0.05/length(modulesCoexpr))])

# affectedModules
# [1] "blue"       "lightgreen" "purple"    


fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Modules enriched in DE/DV genes ", affectedModules), fileConnSummary)
writeLines(paste("\n"), fileConnSummary)

close(fileConnSummary)
# ##############################################################################
