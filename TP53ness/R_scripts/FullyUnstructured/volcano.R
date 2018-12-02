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
TP53_genes <- ifelse((results.df$padj < .01 | abs(results.df$log2FoldChange) > 2.5) , T, F)
TP53_genes <- rownames(res_tableOE_ordered)[which(TP53_genes)]


res_tableOE_ordered <- results.df[order(results.df$padj), ] 
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