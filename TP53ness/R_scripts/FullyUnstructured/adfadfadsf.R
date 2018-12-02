feature.df %>% arrange(desc(Nutlin.3a)) -> feature.df
gdf <- feature.df[,c(genes,"TP53", "Nutlin.3a")]
gdf1 <- gdf[,which(colnames(gdf) %in% c(names(g), "Nutlin.3a", "TP53"))]
gdf1[,which(colnames(gdf1) %in% c(names(g)))] <- scale(gdf1[, which(colnames(gdf1) %in% c(names(g)))])
gdf1 %>% arrange(desc(Nutlin.3a)) -> gdf1
gene.names <- g[which(!is.na(g))]
colnames(gdf1) <- gene.names[colnames(gdf1)]
gdf1 <- gdf1[,which(!is.na(colnames(gdf1)))]
gdf1$TP53 <- gdf$TP53
gdf1$Nutlin.3a <- gdf$Nutlin.3a
gdf1 <- gdf1[which(gdf1$TP53 == 1),]
heat.map <- Heatmap(
  gdf1[,-c((dim(gdf1)[2] - 1), dim(gdf1)[2])],
  name = "genes",
  km = 1,
  show_row_names = F,
  show_column_names = T,
  cluster_rows = F
  ) +
  Heatmap(
  gdf1$Nutlin.3a,
  name = "Nutlin",
  width = unit(5, "mm"),
  cluster_rows = T
  )  



+
  Heatmap(
    gdf1$TP53,
    name = "TP53",
    col = c("red", "blue"),
    cluster_rows = FALSE
  )
