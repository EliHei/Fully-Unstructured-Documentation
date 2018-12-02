library(dplyr)
library(purrr)
library(reshape)
library(reshape2)
library(pheatmap)
library(Questools)
library(jyluMisc)

################### data 
load("../TP53ness_180716_/data/metadata/patmeta_180504.RData")
load("../TP53ness_180716_/data/RNASeq/ddsrna_180301.RData")
load("../TP53ness_180716_/data/Metabolism/Seahorse_20170822_combat.RData")
load("../TP53ness_180716_/data/FISH/fish_20170604.RData")
screen.files <- list.files("../TP53ness_180716_/data/screen/", pattern = ".*\\.rds")
screens <- screen.files %>% map(function(f) readRDS(paste0("../TP53ness_180716_/data/screen/", f)))
names(screens) <- unlist(screen.files %>% map(function(x) gsub(x, pattern = "\\.rds", replacement = ""))) 
wes <- read.csv("../TP53ness_180716_/data/WES/2018-04-09_wes_formatted.csv")
ampliSeq <- read.csv("../TP53ness_180716_/data/ampliSeq/ampliSeq_180509.csv")
ighv <- read.csv("../TP53ness_180716_/data/IGHV/2018-03-05_IGHV_combined.csv")
wgs <- read.csv("../TP53ness_180716_/data/WGS/2018-04-13_wgs_formatted.csv")
load("../TP53ness_180716_/data/Methylome/2017-07-30_meth435k.RData")


################### preprocess
patTP53 <- patMeta[which(!is.na(patMeta$TP53)),]
rownames(patTP53) <- patTP53$Patient.ID 
patTP53 <- data.preproc(patTP53, levels = 4)
write.csv(patMeta, "../TP53ness_180716_/data/Preprocessed/patTP53.csv")


################### RNASeq
dds <- estimateSizeFactors(dds)
dds <- dds[apply(counts(dds), 1, function(x) any(x > 10)),]
dds.vst <- varianceStabilizingTransformation(dds, blind = TRUE)
exprMat <- assay(dds.vst)
sds <- rowSds(exprMat)
names(sds) <- rownames(exprMat)
exprMat <- exprMat[names(sort(sds, decreasing = TRUE)[1:5000]),]
exprMat <- data.frame(exprMat)
rna.data <- t(exprMat)
write.csv(exprMat, file = "../TP53ness_180716_/data/Preprocessed/exprMat.csv")

################### Seahorse
sea.data <- sea$Viability
sea.data <- data.frame(sea = sea.data, row.names = sea$patientID)
sea.data$sea <- impute.continuous(sea.data$sea)
write.csv(exprMat, file = "../TP53ness_180716_/data/Preprocessed/seahorse.csv")
################### FISH
## very sparse!

################### CPS1000
cps.data <- screens$CPS1000
cps.data <- cps.data %>% group_by(patientID, Drug) %>% summarise(mean(normVal, na.rm = T))
colnames(cps.data)[3] <- "normVal"
cps.data <- dcast(cps.data, patientID ~ Drug)
rownames(cps.data) <- cps.data$patientID
cps.data <- cps.data[,-1]
colnames(cps.data) <- paste0(colnames(cps.data), ".")
write.csv(cps.data, "../TP53ness_180716_/data/Preprocessed/CPS1000.csv")

################### 



################### Binding
inr <- Reduce(intersect, list(rownames(cps.data), rownames(rna.data), patTP53$Patient.ID))
d <- cbind(rna.data[inr,], cps.data[inr,])
rownames(d) <- paste0(rownames(d), (patTP53[inr, "TP53"] == 1))
l <- ggm(data.frame(t(d)), rho = 0.3)
g <- l$graph

data1 <- toVisNetworkData(as(g, "igraph"))
group <- as.data.frame(patTP53[inr, "TP53"])$TP53
names(group) <- inr
data1$nodes$group <- group
visNetwork(nodes = data1$nodes, edges = data1$edges)  %>%
  visOptions(highlightNearest = list(
    enabled = TRUE,
    degree = 1,
    hover = TRUE
  ))




kl.div <- div(data = d, g1 = which(grepl(".*TRUE", rownames(d))), g2= which(grepl(".*FALSE", rownames(d))))




library(org.Hs.eg.db)
cols <- c("SYMBOL", "GENENAME")
ensids <- rownames(kl.div)[grep("ENS", rownames(kl.div))]
mapping <- select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")
mapping %>% group_by(ENSEMBL) %>% summarize(sample(SYMBOL, 1)) -> mapping
kl.div$feature <- rownames(kl.div)
colnames(mapping)[2] <- "SYMBOL"
kl.div[mapping$ENSEMBL, "SYMBOL"] <- mapping$SYMBOL
kl.div %>% arrange(desc(KL)) -> kl.div
genes <- kl.div$SYMBOL[!is.na(kl.div$SYMBOL)]
drugs <- kl.div$feature[is.na(kl.div$SYMBOL)]
drugs <- drugs[-grep(pattern = "ENS*", x = drugs)]
top.genes <- head(genes, 100)
top.drugs <- head(drugs, 50)
top <- kl.div$feature[1:100]
top <- union(top, top.drugs)
colnames(d) <- kl.div$feature
d.top <- d[, colnames(d) %in% top]
l <- ggm(data.frame(d.top), rho = 0.4)
g <- l$graph

data1 <- toVisNetworkData(as(g, "igraph"))
group <- as.data.frame(patTP53[as_ids(V(g)), "TP53"])$TP53

data1$nodes$group <- group
data1$nodes %>% dplyr::select(id, label, group) -> data1$nodes
visNetwork(nodes = data1$nodes, edges = data1$edges)  %>%
  visOptions(highlightNearest = list(
    enabled = TRUE,
    degree = 1,
    hover = TRUE
  ))

ggm(data.frame(t(d.top)), rho = 0.35)

cutter <- function(vec, levels){
  to.ret <- cut(vec,breaks = seq((min(vec) - 1), (max(vec) + 1), (max(vec) - min(vec) + 2) /
                         levels),
      labels = 1:levels)
  factor(to.ret)
}



library(SuperLearner)
sl_lasso = SuperLearner(Y = , X = d, family = binomial(),
                        SL.library = "SL.glmnet")

Y = as.numeric(unlist(data.frame(patTP53[rownames(d), "TP53"])))


tp35 <- d.top[(patTP53[inr, "TP53"] == 1),]
ggm(data.frame(t(tp35)), rho = 0.5)

rownames(d.top) <- gsub(pattern = "TRUE", replacement = "", rownames(d.top))
rownames(d.top) <- gsub(pattern = "FALSE", replacement = "", rownames(d.top))
ll <- apply(d.top , 2, cutter, levels = 5)
ll<-data.frame(ll, row.names = row.names(d.top))
l <- min.forest(data.frame(ll))
g <- l$graph

############################
prec <- data.preproc(patTP53)
min.forest(prec)
min.forest(t(d))






###########################

lasso.mod <- glmnet(x = data.matrix(d), y= unlist(Y), alpha = .5, family = "binomial")
coef  <- predict(lasso.mod, type = 'coefficients', s = 0.01)
matrix(coef) -> coef
rownames(coef)[-1] <- colnames(d)
coef <- data.frame(coef)
coef$coef <- abs(coef$coef)
tops <- rownames(coef)[which(coef$coef > 0)]
d.top <- d[, colnames(d) %in% tops]
l <- ggm(data.frame(d.top), rho = 0.4)
g <- l$graph
data1 <- toVisNetworkData(as(g, "igraph"))
group <- as.data.frame(patTP53[as_ids(V(g)), "TP53"])$TP53
data1$nodes$group <- group
data1$nodes %>% dplyr::select(id, label, group) -> data1$nodes
visNetwork(nodes = data1$nodes, edges = data1$edges)  %>%
  visOptions(highlightNearest = list(
    enabled = TRUE,
    degree = 1,
    hover = TRUE
  ))

d.cor <- cor(d.top)
pheatmap::pheatmap(d.cor)



#########################

fc <- cluster_louvain(as.undirected(g))
m<- membership(fc)
com <- names(m)[which(m==44)]
d.com <- d[, com]




acc.calc <- function(train_data,
                     test_data,
                     train_data_class,
                     test_data_class) {
  to.ret <- list()
  5:16 %>% map(function(x)
    class::knn(
      train = train_data ,
      test = test_data ,
      cl = train_data_class ,
      k = x
    )) -> to.ret$pred
  to.ret$pred %>% map(function(x)
    sum(x == test_data_class) / length(test_data_class)) -> to.ret$acc
  print(to.ret$acc)
  to.ret$pred %>% map(function(x)
    table(x, test_data_class)) -> to.ret$table
  to.ret$hm <- to.ret$table %>% map(function(x) pheatmap(x))
  to.ret
}

1:10 %>% map(function(x) sample(x = 1:nrow(d.com), size = nrow(d.com) * 0.8 , replace = T)) %>% map(function(x) acc.calc(train_data = d.com[x , ], test_data = d.com[-x , ], train_data_class = Y[x], test_data_class = Y[-x])) -> cv

cv %>% map(function(x) x$acc) -> a
cv %>% map(function(x) x$hm) -> hm
dd  <-  as.data.frame(matrix(unlist(a), nrow=length(unlist(a[1]))))
rownames(dd) <- 5:16
rowMeans(dd)
hm


#####################
library(caret)
library(MASS)
index <- createDataPartition(Y, p=0.70, list=FALSE)
set.seed(1000)
model.lda<-train(x = d.top[index,],y = factor(Y[index]), method = 'spls',metric = "Accuracy")
pred_test<-predict(object = model.lda,newdata = d.top[-index,])
confusionMatrix(pred_test,factor(Y[-index]))



####################
data.preproc(patTP53) -> p1
rownames(p1) <- rownames(patTP53)
l <- min.forest(data.frame(t(p1)))
g <- l$graph



####################
kl.div <- div(data = p1[, -110], g1 = which(patTP53$TP53 == 0), g2 = which(patTP53$TP53 == 1), permute = 1000)





######################
colnames(d) <- kl.div$feature
d1 <- cbind(d, as.numeric(unlist(patTP53[rownames(d), "TP53"])))
d1 <- cbind(d1, p1[rownames(d),])
rownames(d1) <- rownames(d)
d1 <- d1[, -c(2, 5176)]
kl.div2.1 <- div2(d1, d$Nutlin.3a, as.numeric(unlist(patTP53[rownames(d), "TP53"])), permute = 0)
kl.div2.1 %>% arrange(desc(KL)) -> kl.div2.1

features<- colnames(d1)[order(desc(kl.div2.1$KL))]
d1.sel <- d1[,features[1:10]]
d1.sel <- cbind(d1.sel, as.numeric(unlist(patTP53[rownames(d), "TP53"])))
l <- min.forest(data.frame(d1.sel))


for(i in 1:length(features)){

  
}





divi(feature.df, which(feature.df$Nutlin.3a >0.9 & feature.df$TP53 == 1), which(feature.df$Nutlin.3a <0.9 & feature.df$TP53 == 1)) -> kl.divvv
genes <- rownames(kl.divvv)[order(desc(kl.divvv$KL))]
write(genes, "../data/featuers.txt")







edo2 <- gseNCG(geneList, nPerm=10000, pvalueCutoff = 10)
edo <- enrichDGN(ids$ENTREZID)
DOSE::barplot(edo, showCategory=20)


geneList <- df$V1
names(geneList) <- rownames(df)


library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)


p2 <- enrichplot::dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")








kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 10,
               verbose      = FALSE)





pheatmap(feature.df)


kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 10,
               verbose      = FALSE)
dotplot(kk2)





