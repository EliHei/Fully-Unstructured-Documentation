

```{r TP53, message=FALSE, warning=FALSE, paged.print=FALSE}
design(dds) <- ~ TP53

ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, contrast = c("TP53", 1, 2))
results.df <- as.data.frame(results.cell)
rdf <- results.df[results.df$padj <= 0.0001 & abs(results.df$log2FoldChange) > 3, ]
rdf <- rdf[!is.na(rdf$padj),]
```

```{r TP53, message=FALSE, warning=FALSE, paged.print=FALSE}
design(dds) <- ~ Methylation_Cluster + TP53
ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, contrast = c("TP53", 1, 2))
results.df <- as.data.frame(results.cell)
rdf <- results.df[results.df$padj <= 0.0001 & abs(results.df$log2FoldChange) > 3, ]
rdf <- rdf[!is.na(rdf$padj),]
```

```{r bar_plot}
res <- createInput(results.df, pCut = 0.1, ifFDR = TRUE)
res <- ens2sym(res)
plots <- names(gmts) %>% map(function(x) barplot(res, gmts$x, x))
plots
```


```{r message=FALSE, warning=FALSE, paged.print=FALSE}

design(dds) <- ~ nut.res + Methylation_Cluster
ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, contrast = c("nut.res", 1, 2))
results.df <- as.data.frame(results.cell)
rdf <- results.df[results.df$padj <= 0.0001 & abs(results.df$log2FoldChange) > 3, ]
rdf <- rdf[!is.na(rdf$padj),]
```

```{r message=FALSE, warning=FALSE, paged.print=FALSE}
design(dds) <- ~ nut + Methylation_Cluster
ddsCell <- DESeq(dds)
results.cell <- results(ddsCell, contrast = c("nut", 1, 2))
results.df <- as.data.frame(results.cell)
rdf <- results.df[results.df$padj <= 0.001 & abs(results.df$log2FoldChange) > 3, ]
rdf <- rdf[!is.na(rdf$padj),]
```


```{r}

```



