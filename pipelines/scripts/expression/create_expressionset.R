suppressPackageStartupMessages(library(Biobase))

args <- commandArgs(trailingOnly = TRUE)
eFile <- args[1]
fFile <- args[2]
pFile <- args[3]

eMatrix <- as.matrix(read.table(eFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
fData <- read.table(fFile, row.names=1, header=TRUE, sep="\t")
pData <- read.table(pFile, row.names=1, header=TRUE, sep="\t")

featureData <- new("AnnotatedDataFrame", data=fData)
phenoData <- new("AnnotatedDataFrame", data=pData)
exampleSet <- ExpressionSet(assayData=eMatrix, phenoData=phenoData, featureData=featureData)

saveRDS(exampleSet, "expression_set.rds")