library(org.Hs.eg.db)
library(cogena)
library(tidyverse)
library(Questools)
library(jyluMisc)

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


g1 <- which(feature.df$Nutlin.3a < 0.9)
g2 <- which(feature.df$Nutlin.3a > 0.9)
divi(feature.df, g1, g2) -> kl.divvv
sig.list <- 1:dim(feature.df)[2] %>% map(function(i) sign(mean(feature.df[g1, i]) - mean(feature.df[g2, i])))
sigs <- unlist(sig.list)
kl.divvv$KL <- kl.divvv$KL * sigs
kl.div <- kl.divvv
cols <- c("SYMBOL", "GENENAME")
ensids <- rownames(kl.div)[grep("ENS", rownames(kl.div))]
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")
mapping %>% group_by(ENSEMBL) %>% summarize(sample(SYMBOL, 1)) -> mapping
kl.div$feature <- rownames(kl.div)
colnames(mapping)[2] <- "SYMBOL"
kl.div[mapping$ENSEMBL, "SYMBOL"] <- mapping$SYMBOL
kl.div %>% arrange(desc(KL)) -> kl.div
kl.div %>% filter(!is.na(kl.div$SYMBOL)) -> kl.Symbol
kl.Symbol %>% group_by(SYMBOL) %>% summarise(mean(KL)) -> kls
rownames(kls) <- kls$SYMBOL
my.symbols <- kls$SYMBOL
hs <- org.Hs.eg.db
AnnotationDbi::select(hs, 
                      keys = my.symbols,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "SYMBOL") -> ids
ids %>% group_by(SYMBOL) %>% summarise(unique(ENTREZID)[1]) -> ids
kl <- kls$`mean(KL)`
names(kl) <- kls$SYMBOL
df <- data.frame(V1=sort(kl, decreasing=TRUE))

pdf("results_bar_nutlin.pdf")
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




