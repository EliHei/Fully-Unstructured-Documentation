div <- function(data,
                g1,
                g2,
                permute = 0,
                levels = 5) {
  is.cat <- function(var) {
    !length(unique(var)) > levels
  }
  kl.calc <- function(data, g1, g2) {
    1:dim(data)[2] %>% map(function(x)
      freq(data[, x], g1, g2))  %>% map(function(x)
        abs(KL.plugin(x$g1, x$g2)) + abs(KL.plugin(x$g2, x$g1))) -> to.ret
    unlist(to.ret)
  }
  freq <- function(vec, g1, g2) {
    if (!is.cat(vec))
      vec <-
        cut(vec,
            breaks = seq((min(vec) - .0000001), (max(vec) + .0000001), (max(vec) - min(vec) + .0000002) /
                           levels),
            labels = 1:levels)
    to.ret <- list(g1 = c(), g2 = c())
    levels(factor(vec)) %>% map(function(x)
      list(g1 = max(1, sum(vec[g1] == x)), g2 = max(1, sum(vec[g2] == x)))) %>% map(function(x)
        to.ret <<-
          list(
            g1 = c(to.ret$g1, x$g1),
            g2 = c(to.ret$g2, x$g2)
          )) -> na
    to.ret
  }
  p.val <- function(x, vec) {
    which(sort(vec, decreasing = T) < x)[1] / length(vec)
  }
  data <- data.frame(data)
  g1.g2 <- c(g1, g2)
  kl <-
    kl.calc(data, g1.g2[1:length(g1)], g1.g2[(length(g1) + 1):length(g1.g2)])
  if (permute > 0) {
    kl.df <- data.frame()
    1:permute %>% map(function(x)
      shuffle(g1.g2)) %>% map(function(x)
        list(g1 = x[1:length(g1)], g2 = x[(length(g1) + 1):length(x)])) %>% map(function(f)
          kl.calc(data, f[[1]], f[[2]])) %>% map(function(x)
            kl.df <<- rbind(kl.df, x)) -> na
    
    1:dim(kl.df)[2] %>% map(function(i)
      p.val(kl[i], kl.df[, i])) -> kls
    return(data.frame(
      KL = kl,
      row.names = colnames(data),
      p.value = unlist(kls)
    ))
  }
  return(data.frame(KL = kl, row.names = colnames(data)))
}
