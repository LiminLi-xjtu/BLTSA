mylabelling_fix <- function(X, root = NA, k = 50){

  ndim = 1
  #Compute neighborhood indices
  print('Find nearest neighbors...')
  n <- nrow(X)
  d <- ncol(X)
  D <- DISTMAT(X)
  ind <- matrix(nr=nrow(D),nc=ncol(D))
  for (i in 1:nrow(D)){
    ind[i,] <- order(D[i,])
  }
  ni <- ind[,1:(k+1)]
  I <- array(list(), n)
  for (i in 1:n){
    Ii = ni[i,]
    Ii = Ii[which(Ii!=0)]
    I[[i]]=Ii
  }
  
  
  sigma <- matrix(nr=nrow(X),nc=1)
  U <- array(list(), nrow(X))
  condir = rep(NA,times=n)
  
  
  for (i in 1:n) {
    # Compute correlation matrix W
    Ii = I[[i]]
    kt = length(Ii)
    Xi = (scale(X[Ii, ], center = TRUE, scale = FALSE))
    #Compute local information by computing d largest eigenvectors of W
    svd = svd(t(Xi))
    Ui<-svd[["u"]]
    U[[i]]<-Ui[,1:ndim]
    si = svd[["d"]]
    sigma[i,] = si[length(si)]^2/sum(si^2)
    
    
    Z = si %*% t(svd[["v"]])
    Zii = Z - rep(Z[1],times=kt)
    W = t(Zii)  %*% Zii
    condir[i] = length(W[W>0])/((kt)*(kt-1))
    
  }
  
  End = which(condir>0.95)
  
  
  #determine nonbranch points and label them
  cutoff = quantile(sigma,0.8)
  NonBranch0 = which(sigma<cutoff)
  Branch0 =  setdiff(1:n,NonBranch0)
  NonBranch = NonBranch0
  Branch = Branch0
  
  End = setdiff(End,Branch)
  
  if (is.na(root)){
    root = End[1]
  }
  
  
  
  DD0 = get.knn(X, 40, algorithm=c("cover_tree"))
  
  
  DD0index <- DD0[["nn.index"]]
  DD0dist <- DD0[["nn.dist"]]
  p<-c(rep(1:nrow(DD0index),40))-1 
  q<-c(DD0index)-1 
  x<-c(DD0dist)
  n1<-nrow(DD0index)
  M<- new("dgTMatrix",
          i = as.integer(p),
          j = as.integer(q), x=x, Dim=as.integer(c(n1,n1)))
  M <- as.matrix(M)
  M <- (M+t(M))/2
  DD0<-as(M, "dgTMatrix")
  cidx0 = components(graph.adjacency(DD0[End,End], mode="undirected", weighted=TRUE))
  cidx0 = cidx0[["membership"]]
  num_branch = length(unique(cidx0))
  
  
  DD = get.knn(X[NonBranch,], 10, algorithm=c("cover_tree"))
  DDindex <- DD[["nn.index"]]
  DDdist <- DD[["nn.dist"]]
  p<-c(rep(1:nrow(DDindex),10))-1 
  q<-c(DDindex)-1 
  x<-c(DDdist)
  n1<-nrow(DDindex)
  DD<- new("dgTMatrix",
           i = as.integer(p),
           j = as.integer(q), x=x, Dim=as.integer(c(n1,n1)))
  
  
  cidx <- sc(DD, median(c(as.matrix(DD))), num_branch)  
  
  
  
  label <- matrix(0,nr=nrow(X),nc=1)
  label[NonBranch,]=cidx
  
  S<-array(list(), num_branch)
  for (t in 1:num_branch){
    st=which(label==t)
    S[[t]]=st
  }
  
  #extend the labelling to branch points
  dt <- sigma[Branch,]
  Ib <- order(dt)
  Branch_sort <- Branch[Ib]
  
  
  for (i in 1:length(Branch_sort)){
    
    # step 1: find the easiest branch point
    branch = Branch_sort[i]
    
    # step2: classify current_branch
    
    dist1<-rep(NA,num_branch)
    dist2<-rep(NA,num_branch)
    for (t in 1:num_branch){
      if(length(na.omit(S[[t]])) > 4){
        num = 5
      }else{
        num = length(na.omit(S[[t]]))
      }
      
      dt=D[branch,S[[t]]]
      dt <- na.omit(dt)
      Dist1=sort(dt)
      Ib=order(dt)
      Dist1=Dist1[1:num]
      Ib=Ib[1:num]
      dist1[t]=mean(Dist1)
      near = S[[t]][Ib]
      
      Dist2<-rep(NA,num)
      for (j in 1:num){
        T = U[[near[j]]]
        Dist2[j] = norm((diag(d)-T %*% t(T)) %*% (t(X[branch,])-t(X[near[j],])),"2")
      }
      dist2[t] = mean(Dist2)
      
      
      
      
    }
    
    
    
    
    dist = 0.5*dist1+0.5*dist2
    
    bb <- which.min(c(dist))
    label[branch,] <- bb
    
    e1 <- c(na.omit(S[[bb]]),branch)
    e <- rep(NA,nrow(X))
    e[1:min(length(e1),nrow(X))] <- e1[1:min(length(e1),nrow(X))]
    S[[bb]] <- e
    
    # update branch and nonbranch 
    Branch <- setdiff(Branch_sort, branch)
    NonBranch <- c(NonBranch,branch)
    
    # update tagent space
    Ii <- c(intersect(I[[branch]],S[[bb]]),rep(NA,length(I[[branch]])-length(intersect(I[[branch]],S[[bb]]))))
    Ii <- na.omit(Ii)
    kt <- length(Ii)
    if(length(Ii) == 1){
      Xi <- Xi
    }else{
      Xi <- X[Ii,] - t(matrix(data=colMeans(X[Ii,]), nrow = length(colMeans(X[Ii,])), ncol = length(Ii)))
    }
    
    # Compute local information by computing d largest eigenvectors of W
    svd = svd(t(Xi))
    si = svd[["d"]]
    Ui<-svd[["u"]]
    U[[branch]]<-Ui[,1:ndim]  
    
  }
  
  
  
  result = list()
  result$label = label
  result$cutoff = cutoff
  result$root = root 
  result$condir = condir
  result$sigma = sigma
  return(result)
  
  
}




