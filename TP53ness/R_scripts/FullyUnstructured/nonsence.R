featuers <- colnames(feature.df)
genes <- featuers[grepl("ENS", featuers)]
gdf <- feature.df[,c(genes,"TP53", "Nutlin.3a")]

most.imp <- function(gdf, cor.co = 0.05, t.co = 0.05){
  TP53 <- gdf$TP53
  nut <- gdf$Nutlin.3a
  genes <- colnames(gdf)[-c((dim(gdf)[2]-1), dim(gdf)[2])]
  res <- genes %>% map(function(x) (t.test(gdf[which(TP53 == 1), x], gdf[which(TP53 == 2), x])$p.value < t.co & cor.test(gdf[which(TP53 == 1),x], nut[which(TP53 == 1)])$p.value < cor.co))      
  to.ret <- genes[unlist(res)]
  ens2sym(to.ret)
}

g <- most.imp(gdf, cor.co = 0.05, t.co = 0.01)

