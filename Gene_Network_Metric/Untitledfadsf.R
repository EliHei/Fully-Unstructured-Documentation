library(doParallel)
registerDoParallel(2)
registerDoMC(3)



dist.matrix <- function(df, m, k){
  trans <- matrix.power(t.matrix, k)
  dist.pair <- function(vec1, vec2){
    sqrt(sum((vec1 %*% trans - vec2 %*% trans) ^ 2))
  }
  # to.ret <- 1:dim(df)[2] %>% map(function(x) 1:dim(df)[2] %>% map(function(y) dist.pair(df[,x], df[,y])))
  mat <-  matrix(data=NA, nrow=dim(df)[2], ncol=dim(df)[2])
  mat <- foreach(i = 1:dim(df)[2]) %dopar% {
    # print(i)
    for(j in 1:i){
      # print(i)
      mat[i, j] <<- dist.pair(df[,i], df[,j])
    }
    mat
  }
  mat <- mat %>% map(as.vector)
  mat <- Reduce(function(x,y) pmin(x, y, na.rm = TRUE), mat)
  mat <- matrix(mat, nrow = dim(df)[2], ncol = dim(df)[2])
  mat <- reflect(mat)
  
  # print(m)
  # mat <- reflect(m)
  # diag(mat) <- 0
  mat
}


dist.matrix <- function(df, m=1, k){
  trans <- matrix.power(t.matrix, k)
  dist.pair <- function(vec1, vec2){
    sqrt(sum((vec1 %*% trans - vec2 %*% trans) ^ 2))
  }
  # to.ret <- 1:dim(df)[2] %>% map(function(x) 1:dim(df)[2] %>% map(function(y) dist.pair(df[,x], df[,y])))
  # mat <-  matrix(data=NA, nrow=dim(df)[2], ncol=dim(df)[2])
  mat <- c()
  mat <- foreach(i = 1:dim(df)[2]) %dopar% {
    # print(i)
    for(j in 1:i){
      # print(i)
      mat[j] <<- dist.pair(df[,i], df[,j])
    }
    c(mat, rep(NA, dim(df)[2] - i))
  }
  mat <- as.data.frame(do.call(rbind, mat))
  # mat <- mat %>% map(as.vector)
  # mat <- Reduce(function(x,y) pmin(x, y, na.rm = TRUE), mat)
  # mat <- matrix(mat, nrow = dim(df)[2], ncol = dim(df)[2])
  mat <- reflect(mat)
  # print(m)
  # mat <- reflect(m)
  diag(mat) <- 0
  mat
}

dists <- dist.matrix(hpf[hvg,1:100], k = 5)




trans <- matrix.power(t.matrix, 5)
transf1 <- hpf[hvg,1] %*% trans
transf1 <- unname(transf1)
transf1 <- as.vector(transf1)
transf2 <- hpf[hvg,2] %*% trans
transf2 <- unname(transf2)
transf2 <- as.vector(transf2)
