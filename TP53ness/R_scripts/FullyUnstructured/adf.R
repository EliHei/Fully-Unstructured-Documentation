box.plot <- function(df, x) {
  ggplot(df, aes(x = TP53, y = AGAP4), fill = TP53, col = TP53) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(x) +
    geom_point(aes(col = TP53), size = 10, position = position_jitterdodge()) +
    theme_grey(base_size = 30) 
}

box.plot(data.frame("AGAP4"= feature.df$ENSG000000188234, "del17p" = factor(feature.df[,"del17p"])), "Nutlin | TP53")
fdf <- feature.df[which(feature.df$TP53 == 1),]
df1 <- data.frame(nut = fdf$Nutlin.3a, SCPEP1 = fdf$ENSG00000121064)

Questools::plot(df1, vars = c("nut", "SCPEP1"), levels = 1)
plot(df1,)

cor.test(df1[which(df1$nut < 0.9),]$SCPEP1,  df1[which(df1$nut > 0.9),]$SCPEP1) 


t.test(feature.df[which(feature.df$TP53 == 1),"ENSG00000240563"], feature.df[which(feature.df$TP53 == 2),"ENSG00000240563"])



