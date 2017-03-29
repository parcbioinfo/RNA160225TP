library(foreach)
library(doMC)
registerDoMC()
library(proxy)
library(WGCNA)
library(biomaRt)
library(GOstats)
library("org.Mm.eg.db")
library("edgeR")
library("ecodist")
library("boot")
library("BSDA")

#source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
#biocLite("GOstats")
disableWGCNAThreads()

registerDoMC()

getDoParWorkers()
options(cores=4)
getDoParWorkers()


setwd("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/GeneNetworkAnalysis")

source("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/functionDefinitions.R")

load("../selectedData.RData")

# Get sample id, gender, and group and better format the data frame
sample_info=read.csv("../../samples.txt",header=F)

# divide the data in different groups
HDID_M=as.matrix(unlist(selectedGeneCounts[,sample_info["group",]=="HDID2" & sample_info["gender",]=="M"]))
HDID_F=selectedGeneCounts[,sample_info["group",]=="HDID2" & sample_info["gender",]=="F"]
HSNPT_M=selectedGeneCounts[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="M"]
HSNPT_F=selectedGeneCounts[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="F"]

load("data/adjCoexprModules.RData")
modules=names(table(colorsCoexpr))
########################################################################################################
adj_HDID_M[is.na(adj_HDID_M)]=0
adj_HDID_F[is.na(adj_HDID_F)]=0
adj_HSNPT_M[is.na(adj_HSNPT_M)]=0
adj_HSNPT_F[is.na(adj_HSNPT_F)]=0


######################################################################################################3
fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

HDID_M_vs_HDID_F=c()

wholeNetMantel=mantelCI(1-adj_HDID_M, 1-adj_HDID_F, nboot=30, ncores=10)
wholeNetMantelString=paste(c(as.character(signif(wholeNetMantel, 3))), collapse=" ")
writeLines(paste("\nWhole network Mantel HDID_M vs HDID_F ", wholeNetMantelString, sep=' ', collapse=" "), fileConnSummary)

HDID_M_vs_HDID_F=rbind(HDID_M_vs_HDID_F, wholeNetMantel)

print(paste("\nWhole network Mantel HDID_M vs HDID_F ", wholeNetMantelString, sep=' ', collapse=" "))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleMantel=mantelCI(1-adj_HDID_M[moduleGenes, moduleGenes], 1-adj_HDID_F[moduleGenes, moduleGenes], nboot=30, ncores=10)
  moduleMantelString=paste(c(as.character(signif(moduleMantel, 3))), collapse=" ")
  writeLines(paste("\n", module, " Mantel HDID_M vs HDID_F ", moduleMantelString, sep=' ', collapse=" "), fileConnSummary)
  print(paste("\n", module, " Mantel HDID_M vs HDID_F ", moduleMantelString, sep=' ', collapse=" "))
  HDID_M_vs_HDID_F=rbind(HDID_M_vs_HDID_F, moduleMantel)
  
}
###################################################################################3
wholeNetMantel=mantelCI(1-adj_HSNPT_M, 1-adj_HSNPT_F, nboot=30, ncores=10)
wholeNetMantelString=paste(c(as.character(signif(wholeNetMantel, 3))), collapse=" ")

HSNPT_M_vs_HSNPT_F=c()
HSNPT_M_vs_HSNPT_F=rbind(HSNPT_M_vs_HSNPT_F, wholeNetMantel)

writeLines(paste("\nWhole network Mantel HSNPT_M vs HSNPT_F ", wholeNetMantelString, sep=' ', collapse=" "), fileConnSummary)
print(paste("\nWhole network Mantel HSNPT_M vs HSNPT_F ", wholeNetMantelString, sep=' ', collapse=" "))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleMantel=mantelCI(1-adj_HSNPT_M[moduleGenes, moduleGenes], 1-adj_HSNPT_F[moduleGenes, moduleGenes], nboot=30, ncores=10)
  moduleMantelString=paste(c(as.character(signif(moduleMantel, 3))), collapse=" ")
  writeLines(paste("\n", module, " Mantel HSNPT_M vs HSNPT_F ", moduleMantelString, sep=' ', collapse=" "), fileConnSummary)
  HSNPT_M_vs_HSNPT_F=rbind(HSNPT_M_vs_HSNPT_F, moduleMantel)
  
}

###################################################################################3

wholeNetMantel=mantelCI(1-adj_HDID_M, 1-adj_HSNPT_M, nboot=30, ncores=10)
wholeNetMantelString=paste(c(as.character(signif(wholeNetMantel, 3))), collapse=" ")

HDID_M_vs_HSNPT_M=c()
HDID_M_vs_HSNPT_M=rbind(HDID_M_vs_HSNPT_M, wholeNetMantel)

writeLines(paste("\nWhole network Mantel HDID_M vs HSNPT_M ", wholeNetMantelString, sep=' ', collapse=" "), fileConnSummary)
print(paste("\nWhole network Mantel HDID_M vs HSNPT_M ", wholeNetMantelString, sep=' ', collapse=" "))


for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleMantel=mantelCI(1-adj_HDID_M[moduleGenes, moduleGenes], 1-adj_HSNPT_M[moduleGenes, moduleGenes], nboot=30, ncores=10)
  moduleMantelString=paste(c(as.character(signif(moduleMantel, 3))), collapse=" ")
  writeLines(paste("\n", module, " Mantel HDID_M vs HSNPT_M ", moduleMantelString, sep=' ', collapse=" "), fileConnSummary)
  print(paste("\n", module, " Mantel HDID_M vs HSNPT_M ", moduleMantelString, sep=' ', collapse=" "))
  HDID_M_vs_HSNPT_M=rbind(HDID_M_vs_HSNPT_M, moduleMantel)
  
}
###################################################################################3
wholeNetMantel=mantelCI(1-adj_HDID_F, 1-adj_HSNPT_F, nboot=30, ncores=10)
wholeNetMantelString=paste(c(as.character(signif(wholeNetMantel, 3))), collapse=" ")

HDID_F_vs_HSNPT_F=c()
HDID_F_vs_HSNPT_F=rbind(HDID_F_vs_HSNPT_F, wholeNetMantel)


writeLines(paste("\nWhole network Mantel HDID_F vs HSNPT_F ", wholeNetMantelString, sep=' ', collapse=" "), fileConnSummary)
print(paste("\nWhole network Mantel HDID_F vs HSNPT_F ", wholeNetMantelString, sep=' ', collapse=" "))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleMantel=mantelCI(1-adj_HDID_F[moduleGenes, moduleGenes], 1-adj_HSNPT_F[moduleGenes, moduleGenes], nboot=30, ncores=10)
  moduleMantelString=paste(c(as.character(signif(moduleMantel, 3))), collapse=" ")
  writeLines(paste("\n", module, " Mantel HDID_F vs HSNPT_F ", moduleMantelString, sep=' ', collapse=" "), fileConnSummary)
  print(paste("\n", module, " Mantel HDID_F vs HSNPT_F ", moduleMantelString, sep=' ', collapse=" "))
  HDID_F_vs_HSNPT_F=rbind(HDID_F_vs_HSNPT_F, moduleMantel)
  
}
close(fileConnSummary)


combinedResultsMantel=cbind(HDID_M_vs_HDID_F, HSNPT_M_vs_HSNPT_F, HDID_M_vs_HSNPT_M, HDID_F_vs_HSNPT_F)
combinedResultsMantel=signif(combinedResultsMantel, 3)
colnames(combinedResultsMantel)=c("HDID_M_vs_HDID_F_Mantel", "HDID_M_vs_HDID_F_LI" , "HDID_M_vs_HDID_F_UI" ,
                                  "HSNPT_M_vs_HSNPT_F_Mantel", "HSNPT_M_vs_HSNPT_F_LI" , "HSNPT_M_vs_HSNPT_F_UI",
                                  "HDID_M_vs_HSNPT_M_Mantel", "HDID_M_vs_HSNPT_M_LI" , "HDID_M_vs_HSNPT_M_UI",
                                  "HDID_F_vs_HSNPT_F_Mantel", "HDID_F_vs_HSNPT_F_LI" , "HDID_F_vs_HSNPT_F_UI"
)

rownames(combinedResultsMantel)=c("whole network", modules)

write.csv(combinedResultsMantel, file="resultsCoexpr/combinedResultsMantel.csv")

###################################################################################3

conn_HDID_M=intramodularConnectivity(adj_HDID_M, colorsCoexpr, scaleByMax=F)
conn_HDID_F=intramodularConnectivity(adj_HDID_F, colorsCoexpr, scaleByMax=F)

wholeNetTtest=t.test(conn_HDID_M[,"kTotal"], conn_HDID_F[,"kTotal"], paired = T )

HDID_M_vs_HDID_F=c()
HDID_M_vs_HDID_F=rbind(HDID_M_vs_HDID_F, c(wholeNetTtest$statistic, wholeNetTtest$p.value))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleTtest=t.test(conn_HDID_M[moduleGenes,"kWithin"], conn_HDID_F[moduleGenes,"kWithin"], paired = T )
  HDID_M_vs_HDID_F=rbind(HDID_M_vs_HDID_F, c(moduleTtest$statistic, moduleTtest$p.value))
  
}
###################################################################################3

conn_HSNPT_M=intramodularConnectivity(adj_HSNPT_M, colorsCoexpr, scaleByMax=F)
conn_HSNPT_F=intramodularConnectivity(adj_HSNPT_F, colorsCoexpr, scaleByMax=F)

wholeNetTtest=t.test(conn_HSNPT_M[,"kTotal"], conn_HSNPT_F[,"kTotal"], paired = T )

HSNPT_M_vs_HSNPT_F=c()
HSNPT_M_vs_HSNPT_F=rbind(HSNPT_M_vs_HSNPT_F, c(wholeNetTtest$statistic, wholeNetTtest$p.value))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleTtest=t.test(conn_HSNPT_M[moduleGenes,"kWithin"], conn_HSNPT_F[moduleGenes,"kWithin"], paired = T )
  HSNPT_M_vs_HSNPT_F=rbind(HSNPT_M_vs_HSNPT_F, c(moduleTtest$statistic, moduleTtest$p.value))
  
}


###################################################################################3

conn_HDID_M=intramodularConnectivity(adj_HDID_M, colorsCoexpr, scaleByMax=F)
conn_HSNPT_M=intramodularConnectivity(adj_HSNPT_M, colorsCoexpr, scaleByMax=F)

wholeNetTtest=t.test(conn_HDID_M[,"kTotal"], conn_HSNPT_M[,"kTotal"], paired = T )

HDID_M_vs_HSNPT_M=c()
HDID_M_vs_HSNPT_M=rbind(HDID_M_vs_HSNPT_M, c(wholeNetTtest$statistic, wholeNetTtest$p.value))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleTtest=t.test(conn_HDID_M[moduleGenes,"kWithin"], conn_HSNPT_M[moduleGenes,"kWithin"], paired = T )
  HDID_M_vs_HSNPT_M=rbind(HDID_M_vs_HSNPT_M, c(moduleTtest$statistic, moduleTtest$p.value))
  
}

###################################################################################3

conn_HDID_F=intramodularConnectivity(adj_HDID_F, colorsCoexpr, scaleByMax=F)
conn_HSNPT_F=intramodularConnectivity(adj_HSNPT_F, colorsCoexpr, scaleByMax=F)

wholeNetTtest=t.test(conn_HDID_F[,"kTotal"], conn_HSNPT_F[,"kTotal"], paired = T )

HDID_F_vs_HSNPT_F=c()
HDID_F_vs_HSNPT_F=rbind(HDID_F_vs_HSNPT_F, c(wholeNetTtest$statistic, wholeNetTtest$p.value))

for (module in modules){
  print(module)
  moduleGenes=geneNames[colorsCoexpr==module]
  moduleTtest=t.test(conn_HDID_F[moduleGenes,"kWithin"], conn_HSNPT_F[moduleGenes,"kWithin"], paired = T )
  HDID_F_vs_HSNPT_F=rbind(HDID_F_vs_HSNPT_F, c(moduleTtest$statistic, moduleTtest$p.value))
  
}

########################################################################################################3
combinedResultsTtest=cbind(HDID_M_vs_HDID_F, HSNPT_M_vs_HSNPT_F, HDID_M_vs_HSNPT_M, HDID_F_vs_HSNPT_F)
combinedResultsTtest=signif(combinedResultsTtest, 3)
colnames(combinedResultsTtest)=c("HDID_M_vs_HDID_F_Tstatistic", "HDID_M_vs_HDID_F_pValue" ,
                                 "HSNPT_M_vs_HSNPT_F_Tstatistic", "HSNPT_M_vs_HSNPT_F_pValue" , 
                                 "HDID_M_vs_HSNPT_M_Tstatistic", "HDID_M_vs_HSNPT_M_pValue" , 
                                 "HDID_F_vs_HSNPT_F_Tstatistic", "HDID_F_vs_HSNPT_F_pValue" 
)

rownames(combinedResultsTtest)=c("whole network", modules)

write.csv(combinedResultsTtest, file="resultsCoexpr/combinedResultsTtest.csv")





# plot(coexprConn1[,"kTotal"], coexprConn2[,"kTotal"], xlim=c(0,600), ylim=c(0,600))
# plot(coexprConn1[colorsCoexpr=="white","kWithin"], coexprConn2[colorsCoexpr=="white","kWithin"], xlim=c(0,50), ylim=c(0,50))
# plot(coexprConn1[colorsCoexpr=="white","kTotal"], coexprConn2[colorsCoexpr=="white","kTotal"], xlim=c(0,350), ylim=c(0,350))
# 
# res=t.test(coexprConn1[,"kTotal"], coexprConn2[,"kTotal"], paired=T)
# 
# for (module in modules){
#   moduleGenes=geneNames[colorsCoexpr==module]
#   
#   res=t.test(coexprConn1[moduleGenes,"kWithin"], coexprConn2[moduleGenes,"kWithin"], paired=T)
#   print(module)
#   print(res)
# }
# 
# 

#############################################################################
# start.time <- Sys.time()
# 
# diffConnF = evaluateDiffConn(colorsCoexpr, HSNPT_F, HDID_F, softPower=6, nPerm=1000, nCores=10)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# 
# 
# diffConnM = evaluateDiffConn(colorsCoexpr, HSNPT_M, HDID_M, softPower=6, nPerm=1000, nCores=10)
# 
# save(diffConnM, diffConnF, file="data/diffConn.RData")
# 
# load("data/diffConn.RData")
# 
# diffConnM[is.na(diffConnM)]=0
# sum(abs(diffConnM[,1]) > 2)
# 
# 
# diffConnF[is.na(diffConnF)]=0
# sum(abs(diffConnF[,1]) > 2)
# 
# ##############################################################################################
# start.time <- Sys.time()
# 
# diffConnNPT = evaluateDiffConn(colorsCoexpr, HSNPT_F, HSNPT_M, softPower=6, nPerm=1000, nCores=10)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# 
# 
# diffConnHDID = evaluateDiffConn(colorsCoexpr, HDID_F, HDID_M, softPower=6, nPerm=1000, nCores=10)
# 
# save(diffConnNPT, diffConnHDID, file="data/diffConnNPTvsHDID.RData")
# load("data/diffConnNPTvsHDID.RData")
# 
# 
# diffConnHDID[is.na(diffConnHDID)]=0
# sum(abs(diffConnHDID[,1]) > 2)
# 
# 
# diffConnNPT[is.na(diffConnNPT)]=0
# sum(abs(diffConnNPT[,1]) > 2)
# 
# dna.results.dna=test.individual.genes( t(HSNPT_M), t(HDID_M),scores="cor",distance="abs",rescale.scores=F,check.networks = T, num.permutations=100)
# 
# ##############################################################################################
# 
# moduleDisruptionM= moduleDisruption(adj_HDID_M, adj_HSNPT_M, moduleColors, nPerm=1000, nCores=27)
# 
# sum(moduleDisruptionM[,1] > 3)
