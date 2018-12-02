gsc_oncogenic <- loadGSC(gmts[[2]],type = "gmt")

genes <- gsc_oncogenic$gsc$P53_DN.V1_DN
genes <- intersect(genes, top_expr_nut)

pt <- drug_meta[colnames(exprMat),]
pt <- pt[-which(grepl("NA", rownames(pt))),]
pt$annot <- as.numeric(pt$nut) - as.numeric(pt$TP53)
pt <- pt[order(pt$Nutlin.3a),]
pt$annot <- factor(pt$annot)
exprMat_scaled <- scale(t(exprMat))

colAnno <- data.frame(pt[,c("annot", "Nutlin.3a", "del4p")])
rownames(colAnno) <- rownames(pt)
pheatmap::pheatmap(exprMat[top_expr_nut,rownames(pt)], annotation_col = colAnno, cluster_rows=FALSE, cluster_cols=FALSE, labels_col = "")



pheatmap::pheatmap(t(exprMat_scaled)[ML.drug.n.features.cont,rownames(pt)], annotation_col = colAnno, cluster_rows=T, cluster_cols=FALSE, labels_col = "")



pheatmap::pheatmap(t(exprMat_scaled)[intersect(ML.drug.n.features.cont, top_expr_nut),rownames(pt)], annotation_col = colAnno, cluster_rows=T, cluster_cols=FALSE, labels_col = "")

pheatmap::pheatmap(exprMat[ML.IC50,rownames(pt)], annotation_col = colAnno, cluster_rows=FALSE, cluster_cols=FALSE, labels_col = "")


library(progeny)
pathways = progeny(exprMat, scale=FALSE)
gene_expr = getVarianceStabilizedData(dset)
library(airway)
data(airway)
dset = DESeqDataSetFromMatrix(assay(airway),
                              colData=as.data.frame(colData(airway)), design=~dex)
dset = estimateSizeFactors(dset)
dset = estimateDispersions(dset)
gene_expr = getVarianceStabilizedData(dset)






sample.names <- intersect(rownames(drug_meta), colnames(exprMat))
feature.df <- cbind(drug_meta[sample.names,], t(exprMat[,sample.names]))
feature.df <- muvis::data_preproc(feature.df)
feature.df$group <- as.numeric(feature.df$nut) - as.numeric(feature.df$TP53)

# without drugs
fdf <- feature.df[,-which((colnames(feature.df) %in% c("Nutlin.3a", "TP53", "nut", "del17p")))]
fdf <- fdf[,!(colnames(fdf) %in% make.names(colnames(cps.data)))]
fdf <- scale(data.matrix(fdf))
fdf <- fdf[, colSums(is.na(fdf)) != nrow(fdf)]
lambdas <- (seq(0,2000,10)/100000) * 5 + 0.0001
cfs <- 1:200 %>% map(function(x) sel.group(df = fdf, nut = factor(feature.df$group), TP53 = feature.df$TP53, alpha = 0.8, lambda = lambdas[x], spl = 0.75)) %>% map(function(x) data.matrix(x)[,1])
means <- tapply(unlist(cfs), names(unlist(cfs)), mean)
ses <- tapply(unlist(cfs), names(unlist(cfs)), sd)
sel.rates <- tapply(unlist(cfs), names(unlist(cfs)), sel.rate)
ord <- order(means)
dat <- data.frame(feature = names(means)[ord], mean = means[names(means)[ord]], sel.rate = sel.rates[names(means)[ord]], ses = ses[names(means)[ord]])
dat <- dat[dat$mean != 0,]
#dat <- dat[-grep("ENS",dat$feature),]
dat %>% arrange(mean) -> dat
dat$feature <- factor(dat$feature, levels = dat$feature)
ggplot(dat, aes(x=feature, weight=mean, fill = sel.rate/200)) + geom_bar() + coord_flip() +labs(y = "mean coefficient", fill = "selection frequency") 
ML.drug.n.features.cont <- as.character(dat$feature)
write(ML.drug.n.features.cont, "../results/ML_group.txt")
ML.drug.n.features.cont <- read.table("../results/ML_group.txt")
ML.drug.n.features.cont <- ML.drug.n.features.cont[,1]