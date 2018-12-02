library(DESeq2)
library(Questools)
library(dplyr)
library(reshape2)

load("/g/huber/projects/nct/cll/ProcessedData/patAnnotation/patmeta_180504.RData")
load("/g/huber/projects/nct/cll/ProcessedData/RNAseq/objects/ddsrna_180301.RData")
load("/g/huber/projects/nct/cll/ProcessedData/DrugScreens/CPS1000/CPS1000_180503.RData")
load("/g/huber/projects/nct/cll/ProcessedData/Methylome/2017-07-30_meth435k.RData")
load("/g/huber/projects/nct/cll/ProcessedData/Metabolism/Seahorse_20170822.RData")
load("/g/huber/projects/nct/cll/RawData/Methylome/metadata/patmeta.RData")


## Patient Annotaion
patTP53 <- patMeta[which(!is.na(patMeta$TP53)),]
patTP53 <- data.preproc(patTP53, levels = 4)
rownames(patTP53) <- patMeta$Patient.ID[which(!is.na(patMeta$TP53))]
write.csv(patTP53, "data/Preprocessed/patTP53.csv")

## RNASeq
dds <- estimateSizeFactors(dds)
dds <- dds[apply(counts(dds), 1, function(x) any(x > 10)),]
dds.vst <- varianceStabilizingTransformation(dds, blind = TRUE)
exprMat <- assay(dds.vst)
sds <- rowSds(exprMat)
names(sds) <- rownames(exprMat)
exprMat <- exprMat[names(sort(sds, decreasing = TRUE)[1:5000]),]
exprMat <- data.frame(exprMat)
rna.data <- t(exprMat)
write.csv(rna.data, file = "data/Preprocessed/RNAseq.csv")

## Drug
cps.data <- pheno1000
cps.data <- cps.data %>% group_by(patientID, Drug) %>% summarise(mean(normVal, na.rm = T))
colnames(cps.data)[3] <- "normVal"
cps.data <- dcast(cps.data, patientID ~ Drug)
cps.data <- cps.data[-dim(cps.data)[1],]
rownames(cps.data) <- cps.data$patientID
cps.data <- cps.data[,-1]
write.csv(cps.data, "data/Preprocessed/CPS1000.csv")

## Methylome
meth_ids <- read.csv("/g/huber/projects/nct/cll/RawData/Methylome/metadata/2017-06-12_patID_matches.csv")



## Seahorse
sea.data <- sea$Viability
sea.data <- data.frame(sea = sea.data, row.names = sea$patientID)
sea.data <- data.preproc(sea.data)
rownames(sea.data) <- sea$patientID
write.csv(sea.data, file = "data/Preprocessed/seahorse.csv")

## 

