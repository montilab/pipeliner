suppressPackageStartupMessages(library(Biobase))

args <- commandArgs(trailingOnly = TRUE)
eFile <- args[1]
pFile <- args[2]
mFile <- args[3]

exprs <- as.matrix(read.table(eFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
pData <- read.table(pFile, row.names=1, header=TRUE, sep="\t")
mData <- read.table(mFile, row.names=1, header=TRUE, sep="\t")

phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=mData)
exampleSet <- ExpressionSet(assayData=exprs, phenoData=phenoData)

saveRDS(exampleSet, "expression_set.rds")