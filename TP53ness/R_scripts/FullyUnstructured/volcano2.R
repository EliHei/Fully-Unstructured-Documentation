load(file = "../TP53ness_ultimate/results/nut_DEG.RData")
# results.cell <- results(ddsCell, contrast = c("nut.res", "1", "2"))
results.df_nut <- as.data.frame(results.cell)

cols <- c("SYMBOL", "GENENAME")
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(results.df_nut), columns=cols, keytype="ENSEMBL")
mapping %>% dplyr::group_by(ENSEMBL) %>% dplyr::summarize(sample(SYMBOL, 1)) -> mapping
colnames(mapping)[2] <- "SYMBOL"
rownames(mapping) <- mapping$ENSEMBL
results.df_nut$SYMBOL <- mapping[rownames(results.df_nut),]$SYMBOL
colnames(mapping)[2] <- "SYMBOL"
results.df_nut %>% filter(!is.na(SYMBOL)) -> results.df_nut
res %>% dplyr::group_by(SYMBOL) %>% dplyr::summarize(stat = mean(stat)) -> res
results.df_nut <- data.frame(results.df_nut)
results.df_nut <- aggregate(results.df_nut, list(results.df_nut$SYMBOL), mean)
rownames(results.df_nut) <- results.df_nut$Group.1
results.df_nut <- results.df_nut[,-1]
results.df_nut <- results.df_nut[,-7]
nut_genes <- ifelse((results.df_nut$padj < .1 | abs(results.df_nut$log2FoldChange) > 1) , T, F)
nut_genes <- rownames(results.df_nut)[which(nut_genes)]


res_tableOE_ordered <- results.df_nut[order(results.df_nut$padj), ] 
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