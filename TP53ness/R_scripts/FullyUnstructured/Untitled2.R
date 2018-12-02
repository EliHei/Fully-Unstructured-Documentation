library(tidyverse)
library(org.Hs.eg.db)
library(DESeq2)
library("RColorBrewer")
library(ComplexHeatmap)


dds$PatID
dds <- estimateSizeFactors(dds)
dds <- dds[,which(colnames(dds) %in% rownames(feature.df))]
feature.df %>% mutate(nut.res = ifelse(Nutlin.3a > 0.9, 2, 1)) -> feature.df
dds$nut.res <- feature.df$nut.res
design(dds) <- ~ nut.res
ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, contrast = c("nut.res", 1, 2))
results.df <- as.data.frame(results.cell)
rdf <- results.df[results.df$padj <= 0.001 & abs(results.df$log2FoldChange) > 3, ]
rdf <- rdf[!is.na(rdf$padj),]




createInput <- function(DEres, pCut = 0.05, ifFDR = TRUE, rankBy = "stat") {
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
p <- jyluMisc::runGSEA(inputTab = res, gmtFile = gmts$H, GSAmethod = "gsea")

pdf("results_bar_DEG_res.pdf")
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


gmt <- loadGSC(gmts[[1]])
Hypoxia <- gmt$gsc$HALLMARK_HYPOXIA
Hypoxia.genes <- rownames(res)[which(rownames(res) %in% Hypoxia)]
feature.df %>% arrange(desc(Nutlin.3a)) -> feature.df
fens <- feature.df[, colnames(feature.df) %in% mapping$ENSEMBL]
fens <- scale(fens)


fens <- feature.df[,colnames(feature.df) %in% rownames(rdf)]
fens <- scale(fens)
Heatmap(
  fens,
  name = "genes",
  km = 1,
  show_row_names = T,
  show_column_names = T,
  cluster_rows = FALSE
) +
  Heatmap(feature.df$Nutlin.3a,
          name = "Nutlin",
          width = unit(5, "mm"),
          cluster_rows = FALSE) +
  Heatmap(feature.df$Methylation_Cluster,
          name = "Methylation_Cluster",
          col = c("blue", "purple", "red"),
          cluster_rows = FALSE) +
  Heatmap(feature.df$TP53, name = "TP53", col = c("red", "blue"),
          cluster_rows = FALSE) 







