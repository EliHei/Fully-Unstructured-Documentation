library(tidyverse)
library(org.Hs.eg.db)
library(DESeq2)
library("RColorBrewer")
library(ComplexHeatmap)
library(reshape2)
library(gridExtra)
library(DESeq2)
library(Questools)
library(ggrepel)

dds <- estimateSizeFactors(dds)
dds <- dds[,which(dds$diag == "CLL")]
feature.df[which(feature.df$TP53 == 1),] -> fdf
dds <- dds[,which(colnames(dds) %in% rownames(fdf))]
fdf$nut.res <- ifelse(fdf$Nutlin.3a > 0.9 , 2, 1)
dds$nut.res <- fdf[dds$PatID, "nut.res"]
dds$Methylation_Cluster <- fdf[dds$PatID, "Methylation_Cluster"]
design(dds) <- ~ nut.res + Methylation_Cluster
ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, name="nut.res")
save(results.cell, file = "../TP53ness_ultimate/results/nut_DEG.RData")
load(file = "../TP53ness_ultimate/results/nut_DEG.RData")
load(file = "../TP53ness_ultimate/results/TP53_DEG.RData")
# results.cell <- results(ddsCell, contrast = c("nut.res", "1", "2"))
results.df <- as.data.frame(results.cell)

cols <- c("SYMBOL", "GENENAME")
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(results.df), columns=cols, keytype="ENSEMBL")
mapping %>% dplyr::group_by(ENSEMBL) %>% dplyr::summarize(sample(SYMBOL, 1)) -> mapping
colnames(mapping)[2] <- "SYMBOL"
rownames(mapping) <- mapping$ENSEMBL
results.df$SYMBOL <- mapping[rownames(results.df),]$SYMBOL
colnames(mapping)[2] <- "SYMBOL"
results.df %>% filter(!is.na(SYMBOL)) -> results.df
res %>% dplyr::group_by(SYMBOL) %>% dplyr::summarize(stat = mean(stat)) -> res
results.df <- data.frame(results.df)
results.df <- aggregate(results.df, list(results.df$SYMBOL), mean)
rownames(results.df) <- results.df$Group.1
results.df <- results.df[,-1]
results.df <- results.df[,-7]



res_tableOE_ordered <- results.df[order(results.df_nut$padj), ] 
res_tableOE_ordered$genelabels <- ifelse(rownames(res_tableOE_ordered) %in% intersect(TP53_genes, nut_genes) , T, F)
res_tableOE_ordered$threshold <- ifelse((res_tableOE_ordered$padj < .01 | abs(res_tableOE_ordered$log2FoldChange) > 3) , T, F)
# res_tableOE_ordered <- res_tableOE_ordered[-log10(res_tableOE_ordered$padj) < 10,]
ggplot(res_tableOE_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour=threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(res_tableOE_ordered),""), colour=threshold)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + theme_bw()
  
  
  
  
rdf <- results.df[results.df$padj <= 0.1 & abs(results.df$log2FoldChange) > .1, ]
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
mapping %>% dplyr::group_by(ENSEMBL) %>% dplyr::summarize(sample(SYMBOL, 1)) -> mapping
colnames(mapping)[2] <- "SYMBOL"
rownames(mapping) <- mapping$ENSEMBL
res$SYMBOL <- mapping[rownames(res),]$SYMBOL
colnames(mapping)[2] <- "SYMBOL"
res %>% filter(!is.na(SYMBOL)) -> res
res %>% dplyr::group_by(SYMBOL) %>% dplyr::summarize(stat = mean(stat)) -> res
res <- data.frame(res)
rownames(res) <- res$SYMBOL
res %>% dplyr::select(stat) -> res
p <- jyluMisc::runGSEA(inputTab = res, gmtFile = gmts$H, GSAmethod = "gsea")
df <- res
pdf("results_bar_DEG_nutsdf.pdf")
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

gmtfile <- system.file("extdata","h.all.v5.1.symbols.gmt",package="BloodCancerMultiOmics2017")
gmt <- loadGSC(gmtfile)
Hypoxia <- gmt$gsc$HALLMARK_HYPOXIA
Hypoxia.genes <- rownames(res)[which(rownames(res) %in% Hypoxia)]
feature.df %>% arrange(desc(TP53)) -> feature.df
fens <- feature.df[, colnames(feature.df) %in% mapping$ENSEMBL]
fens <- scale(fens)

gmtfile <- system.file("extdata","h.all.v5.1.symbols.gmt",package="BloodCancerMultiOmics2017")
gmt <- loadGSC(gmtfile)
Hypoxia <- gmt$gsc$HALLMARK_HYPOXIA
Hypoxia.genes1 <- rownames(res)[which(rownames(res) %in% Hypoxia)]
feature.df %>% arrange(desc(TP53)) -> feature.df
fens <- feature.df[, colnames(feature.df) %in% mapping$ENSEMBL]
fens <- scale(fens)

require(VennDiagram)
venn.list <- Vennerable::Venn(list(WT_vs_MUT = Hypoxia.genes1[1:33], res_vs_sens = Hypoxia.genes[1:30]))
Vennerable::plot(venn.list, doWeights = F, type = "circles", show = list(Faces = FALSE))


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




ggplot(box, aes(x = factor(year), y = case, fill = code))+
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge())


m.fens <- melt(fens)





genes <- colnames(fens)


cols <- c("SYMBOL", "GENENAME")
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=genes, columns=cols, keytype="ENSEMBL")
mapping %>% group_by(ENSEMBL) %>% summarize(sample(SYMBOL, 1)) -> mapping
colnames(mapping)[2] <- "SYMBOL"
mapping <- data.frame(mapping)
rownames(mapping) <- mapping$ENSEMBL
fens <- fens[, which(colnames(fens) %in% mapping$ENSEMBL)]
colnames(fens) <- mapping[colnames(fens), "SYMBOL"]
genes <- mapping$SYMBOL
genes <- genes[which(!is.na(genes))]

plotter <- function(df, x){
  print(x)
  ggplot(df, aes(x = factor(TP53), y = gene), fill = factor(TP53), col = factor(TP53)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle(x) +
    geom_point(aes(col = factor(TP53)), size = 1, position = position_jitterdodge()) 
}
genes %>% map(function(x) plotter(data.frame("gene"= fens[,x], "TP53" = feature.df$TP53), x)) -> plotList


grid.arrange(grobs=plotList, ncol=2)







ggplot(res_tableOE_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(res_tableOE_ordered),""))) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
