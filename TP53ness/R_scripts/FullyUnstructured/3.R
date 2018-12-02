library(tidyverse)
library(org.Hs.eg.db)
library(DESeq2)

datasets <-
  list.files("/Volumes/elihei/Internship/TP53ness_180716_/data/Preprocessed/") %>% map(function(f)
    read.table(
      paste0("/Volumes/elihei/Internship/TP53ness_180716_/data/Preprocessed/", f),
      sep = ",",
      header = T,
      row.names = 1
    )) 

row.indices <- Reduce(intersect, datasets %>% map(rownames))
feature.df <- do.call(cbind, datasets %>% map(function(x) x[row.indices,]))
feature.df <- feature.df[colSums(!is.na(feature.df)) > 0]
feature.df <- data.preproc(feature.df)
rownames(feature.df) <- row.indices
feature.df <- Filter(function(x)(length(unique(x))>1), feature.df)

dds$PatID
dds <- estimateSizeFactors(dds)
dds <- dds[,which(colnames(dds) %in% rownames(fdf))]

feature.df[which(feature.df$TP53 == 1),] -> fdf
fdf$nut.res <- ifelse(fdf$Nutlin.3a > 0.9 , 2, 1)
dds$nut.res <- fdf$nut.res
design(dds) <- ~ nut.res
ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, contrast = c("nut.res", 1, 2))
results.df <- as.data.frame(results.cell)




createInput <- function(DEres, pCut = 0.05, ifFDR = FALSE, rankBy = "stat") {
  if (ifFDR) {
    inputTab <- data.frame(DEres) %>% rownames_to_column(var = "symbol") %>% filter(padj <= pCut)
  } else {
    inputTab <- data.frame(DEres) %>% rownames_to_column(var = "symbol") %>%filter(pvalue <= pCut)
  }
  
  inputTab <- arrange(inputTab, pvalue) %>% filter(!duplicated(symbol)) %>% select_("symbol", rankBy) %>% data.frame(stringsAsFactors = FALSE)
  rownames(inputTab) <- inputTab$symbol
  inputTab$symbol <- NULL
  colnames(inputTab) <- "stat"
  return(inputTab)
}

gmts <- list(H=system.file("extdata","h.all.v5.1.symbols.gmt",package="BloodCancerMultiOmics2017"),
             C6=system.file("extdata","c6.all.v5.1.symbols.gmt", package="BloodCancerMultiOmics2017"),
             KEGG=system.file("extdata","c2.cp.kegg.v5.1.symbols.gmt", package = "BloodCancerMultiOmics2017"))

res <- createInput(results.df, pCut = 0.1, ifFDR = F)
cols <- c("SYMBOL", "GENENAME")
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(res), columns=cols, keytype="ENSEMBL")
mapping %>% group_by(ENSEMBL) %>% summarize(sample(SYMBOL, 1)) -> mapping
colnames(mapping)[2] <- "SYMBOL"
rownames(mapping) <- mapping$ENSEMBL
res$SYMBOL <- mapping[rownames(res),]$SYMBOL
colnames(mapping)[2] <- "SYMBOL"
res %>% filter(!is.na(SYMBOL)) -> res
res %>% group_by(SYMBOL) %>% summarize(stat = mean(stat)) -> res
res <- data.frame(res)
rownames(res) <- res$SYMBOL
res %>% dplyr::select(stat) -> res

pdf("results_bar_DEG_muVSwt.pdf")
gmtfile <- system.file("extdata","h.all.v5.1.symbols.gmt",package="BloodCancerMultiOmics2017")
rg <- runGSEA(df, gmtfile)
rg$Name <- gsub("HALLMARK_","", rg$Name)
rg <- list(rg)
names(rg) <- "HALLMARK"
enBar <- plotEnrichmentBar(rg, pCut = 0.1, ifFDR = FALSE, setName = "HALLMARK")
plot(enBar)

gmtfile <- system.file("extdata","c6.all.v5.1.symbols.gmt",package="BloodCancerMultiOmics2017")
rg <- runGSEA(df, gmtfile)
rg <- list(rg)
names(rg) <- "ONCOGENIC"
enBar <- plotEnrichmentBar(rg, pCut = 0.1, ifFDR = FALSE, setName = "ONCOGENIC")
plot(enBar)

anno <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
gmtfile <- system.file("extdata", anno, package="cogena")
rg <- runGSEA(df, gmtfile)
rg$Name <- gsub("KEGG_","", rg$Name)
rg <- list(rg)
names(rg) <- "KEGG"
enBar <- plotEnrichmentBar(rg, pCut = 0.1, ifFDR = FALSE, setName = "KEGG")
plot(enBar)
dev.off()


rg


