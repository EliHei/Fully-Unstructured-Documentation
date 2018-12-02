library(ggpubr)
library(limma)
htseq <- read.table("../data/new_batch/GSE81682_HTSeq_counts.txt", header = T)
rownames(htseq) <- htseq$ID
htseq <- htseq[,-1]
cell_type <- read.table("../data/new_batch/all_cell_types.txt", header = T)


cell_type <- cell_type[colnames(htseq),]
cell_type <- cell_type[1:10]
is.one <- apply(cell_type, 1, function(x) names(which(x == 1)))
cell_types <- is.one %>% map(function(x) paste(x, collapse = " "))
cell_types <- unlist(cell_types)
cell_types[which(cell_types == "")] <- "LTHSC_broad"
cell_types[which(grepl("^MPP", cell_types))] <- "MPP"
loop_types <- cell_types[which((cell_types %in% c("CMP_broad", "GMP_broad", "LMPP_broad")) | grepl("^MPP", cell_types))]
loop_types <- names(loop_types)

hvg <- hv.genes(htseq, 300)
g <- ggm1(data = t(htseq[hvg,]), rho = 0.3)
g$graph$network

t.matrix <- abs(g$wi)
t.matrix <- ifelse(t.matrix > 0.05, t.matrix, 0)
diag(t.matrix) <- 0
t.matrix <- data.frame(t.matrix)
outs <- which(colSums(t.matrix, na.rm = T) == 0)
t.matrix <- t.matrix[-outs, -outs]
hvg <- hvg[-outs]
t.matrix <- t(scale(t.matrix, center=FALSE, scale=colSums(t.matrix)))








# d.matrix <- read.csv("../results/eu.csv", row.names = 1)
d.matrix <- dist(t(htseq[hvg,]))
# ot <- colnames(htseq)

dm <- DiffusionMap(data = data.frame(sample = colnames(htseq[,ot]), col = strsplit2(colnames(htseq[,ot]), "_")[,1]), distance = as.dist(d.matrix))
dpt <- DPT(dm)


dm <- DiffusionMap(data = data.frame(sample = colnames(htseq), col = strsplit2(colnames(htseq), "_")[,1]), distance = as.dist(d.matrix))
dpt <- DPT(dm)

# plot.DPT(dpt, dcs = c(1, 2,3), col = factor(meta_data$Batch_desc[CMP]))
df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = dpt$col)
p111<- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr() + 
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = dpt$col)
p112 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = dpt$col)
p113 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()

# d.matrix <- read.csv("../results/eu.csv", row.names = 1)
# ot <- rownames(d.matrix)
# ct <- cell_types[ot]
ct <- cell_types
ct[which(is.na(ct))] <- "LTHSC_broad"
# dm <- DiffusionMap(data = data.frame(sample = colnames(htseq[,ot]), col = ct), distance = as.dist(d.matrix))
# dpt <- DPT(dm)

# plot.DPT(dpt, dcs = c(1, 2,3), col = factor(meta_data$Batch_desc[CMP]))
df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = ct)
p211 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = ct)
p212 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = ct)
p213 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()



d.matrix <- read.csv("../results/ggm_03_new_batch.csv", row.names = 1)

dm <- DiffusionMap(data = data.frame(sample = colnames(htseq), col = strsplit2(colnames(htseq), "_")[,1]), distance = as.dist(d.matrix))
dpt <- DPT(dm)

df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = dpt$col)
p121 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = dpt$col)
p122 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = dpt$col)
p123 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


# d.matrix <- read.csv("../results/ggm_03_new_batch.csv", row.names = 1)
# 
# dm <- DiffusionMap(data = data.frame(sample = colnames(htseq), col = cell_types), distance = as.dist(d.matrix))
# dpt <- DPT(dm)

df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = cell_types)
p221 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = cell_types)
p222 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()

df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = cell_types)
p223 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


d.matrix <- read.csv("../results/euc_htseq.csv", row.names = 1)
ot <- rownames(d.matrix)
d.matrix <- read.csv("../results/scater_htseq_.csv", row.names = 1)

dm <- DiffusionMap(data = data.frame(sample = colnames(htseq[ot]), col = strsplit2(colnames(htseq[ot]), "_")[,1]), distance = as.dist(d.matrix))
dpt <- DPT(dm)

df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = dpt$col)
p131 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = dpt$col)
p132 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = dpt$col)
p133 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()

ct <- cell_types[ot]
ct[which(is.na(ct))] <- "LTHSC_broad"
# dm <- DiffusionMap(data = data.frame(sample = colnames(htseq[,ot]), col = ct), distance = as.dist(d.matrix))
# dpt <- DPT(dm)

# plot.DPT(dpt, dcs = c(1, 2,3), col = factor(meta_data$Batch_desc[CMP]))
df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = ct)
p231 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = ct)
p232 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = ct)
p233 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()



d.matrix <- read.csv("../results/eu.csv", row.names = 1)
ot <- rownames(d.matrix)
d.matrix <- read.csv("../results/hvg_cyto.csv", row.names = 1)

dm <- DiffusionMap(data = data.frame(sample = colnames(htseq[ot]), col = strsplit2(colnames(htseq[ot]), "_")[,1]), distance = as.dist(d.matrix))
dpt <- DPT(dm)

df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = dpt$col)
p141 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = dpt$col)
p142 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = dpt$col)
p143 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


ct <- cell_types[ot]
ct[which(is.na(ct))] <- "LTHSC_broad"
# dm <- DiffusionMap(data = data.frame(sample = colnames(htseq[,ot]), col = ct), distance = as.dist(d.matrix))
# dpt <- DPT(dm)

# plot.DPT(dpt, dcs = c(1, 2,3), col = factor(meta_data$Batch_desc[CMP]))
df <- data.frame(DPT1 = dpt$DC1, DPT2 = dpt$DC2, cell_type = ct)
p241 <- ggplot(df, aes(x = DPT1, y = DPT2, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT2 = dpt$DC2, DPT3 = dpt$DC3, cell_type = ct)
p242 <- ggplot(df, aes(x = DPT2, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


df <- data.frame(DPT1 = dpt$DC1, DPT3 = dpt$DC3, cell_type = ct)
p243 <- ggplot(df, aes(x = DPT1, y = DPT3, color = cell_type)) + geom_point() + theme_pubr()


DPT.fig1 <- ggarrange(p111, p112, p113, p121, p122, p123, p131, p132, p133, p141, p142, p143, nrow= 4, ncol=3,
                     common.legend = TRUE, legend="right", labels = c("Euclidean", "", "", "Reg. GNM", "", "", "Ido GNM", "", "", "EN + Ido GNM", "", "" ))
                     
DPT.fig2 <- ggarrange(p211, p212, p213, p221, p222, p223, p231, p232, p233, p241, p242, p243, nrow= 4, ncol=3, 
                     common.legend = TRUE, legend="right", labels = c("Euclidean", "", "", "Reg. GNM", "", "", "Ido GNM", "", "", "EN + Ido GNM", "", "" ))


DPT.fig2 <- ggarrange(p211, p212, p213, nrow= 1, ncol=3, 
                      common.legend = TRUE, legend="right", labels = c("EN + Euclidean", "", ""))

DPT.fig2 <- ggarrange(p111, p112, p113, nrow= 1, ncol=3, 
                      common.legend = TRUE, legend="right", labels = c("EN + Euclidean", "", ""))



ido <- read.table("../data/Ido_data/GSE72857_umitab.txt")
rownames(ido) <- strsplit2(rownames(ido), ";")[,1]


hpf <- read.table("../data/nestorowa_corrected_log2_transformed_counts.txt", header = T)
htseq <- t(hpf)


hvg <- intersect(hvg, rownames(ido))
outs <- which(rowSds(data.matrix(ido[hvg,])) == 0)
hvg <- hvg[-outs]
g <- ggm1(data = data.frame(t(ido[intersect(hvg, rownames(ido)),])), rho = 0.1)


df_temp <- data.matrix(htseq[hvg, ])
df_temp <- ifelse(df_temp > 0, 1, 0)
df_temp <- data.frame(df_temp)
classes <- strsplit2(colnames(htseq), "_")[,1]
HSPC.df <- df_temp[,which(classes == "HSPC")]
LTHSC.df <- df_temp[,which(classes == "LT.HSC")]
Prog.df <- df_temp[,which(classes == "Prog")]

HSPC_groups <- apply(HSPC.df, 1, mean)
HSPC_size <- apply(HSPC.df, 1, function(x) (1/1+(sd(x)))^4)
LTHSC_groups <- apply(LTHSC.df, 1, mean)
LTHSC_size <- apply(LTHSC.df, 1, function(x) (1/1+(sd(x)))^4)
Prog_groups <- apply(Prog.df, 1, mean)
Prog_size <- apply(Prog.df, 1, function(x) (1/1+(sd(x)))^4)


l <- g$graph$graph
l1 <- delete.isolates(l)
lay <- layout_nicely(l1)
# CMP CD41
graph.col.red(graph = l, lay= lay, grp = (as.numeric(HSPC_groups) + 1), sz = HSPC_size)
title("HSPC",cex.main=1,col.main="black")
# CMP Flt3+ Csf1r+
graph.col.green(graph = l, lay= lay, grp = (as.numeric(LTHSC_groups) + 1), sz = LTHSC_size)
title("LTHSC",cex.main=1,col.main="black")
# CMP Irf8-GFP+ MHCII+
graph.col.blue(graph = l, lay= lay, grp = (as.numeric(Prog_groups) + 1), sz = Prog_size)
title("Progenitor",cex.main=1,col.main="black")


graph.col.green <- function(graph, lay, grp, sz){
  igraph::V(graph)$color <- rgb(0, floor(255*(grp^4)/max((grp^4))), maxColorValue=255, alpha=255, 0)
  v <-  igraph::V(graph)
  graph <- as.undirected(graph)
  igraph::V(graph)$size <- sz
  iso <- function(graph, mode = 'all') {
    isolates <- which(degree(graph, mode = mode) == 0) 
    isolates
  }
  is <- iso(graph)
  print(is)
  graph <- delete.isolates(graph)
  plot.igraph(
    graph,
    vertex.size = V(graph)$size,
    layout = lay[-is,],
    vertex.frame.color = igraph::V(graph)$color
  )
}




V(l)$name = names(V(l))
write.graph(l, "../results/ggm_0.3.gml")
                     
