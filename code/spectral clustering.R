sc <- function(A, sigma, num_clusters){
  
  n <- nrow(A)
  if(sigma == 0){
    col_count <- colSums(as.matrix(A)!=0) 
    col_sum <- colSums(as.matrix(A))
    col_mean <- col_sum/col_count
    
    x <- A@i+1
    y <- A@j+1
    val <- A@x
    val <- -val*val/col_mean[x]/col_mean[y]/2
    A <- new("dgTMatrix",
             i = as.integer(A@i),
             j = as.integer(A@j), x=val, Dim=as.integer(c(n,n)))
  }else{
    A <- A*A
    A <- -A/(2*sigma*sigma)
    
  }
  
  
  num <- 2000
  num_iter <- 300
  S<-matrix(nrow = nrow(A))
  S<-S[,-1]
  for(i in 1:num_iter){
    start_index = 1 + (i-1)*num
    if (start_index>ncol(A)){
      break
    }
    end_index = min(i*num, n)
    S1<-exp(A[,start_index:end_index])
    S1[which(as.matrix(S1) == 1)] <- 0
    S<-cbind(S,S1)
  }
  
  S[is.na(S)] <- 0
  S <- as(S,'dgCMatrix')
  
  D <- rowSums(as.matrix(S))+ (1e-10)
  D <- sqrt(1/D)
  D <- diag(D)
  D <- as(D,'dgCMatrix')
  L <- D %*% S %*% D
  
  eig <- eigs(L, num_clusters,which = "LM")
  V <- eig[["vectors"]]
  
  sq_sum <- sqrt(rowSums(V*V)) + 1e-20
  U <- V/matrix(data=sq_sum, nrow = length(sq_sum), ncol = num_clusters)
  
  cluster_labels <- kmeans(U, num_clusters)
  cluster_labels <- cluster_labels[["cluster"]]
  
  return(cluster_labels)
}
