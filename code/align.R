align <- function(X, label, root, parameters){
  n = nrow(X)
  d = 1
  
  labelset = unique(label)
  labelset = sort(labelset)
  Branch<-array(NA, dim = c(1, nrow(label), length(unique(label))))
  for (t in 1:length(labelset)){
    st=which(label==t)
    if (length(st)>0){
      Branch[1,1:length(st) ,labelset[t]]=st
    }else{
      Branch[1,,labelset[t]]=NA
    }
  }
  
  rootlabel = label[root,]
  nonrootlabel = setdiff(labelset,rootlabel)
  RootBranch = Branch[,,rootlabel]
  RootBranch = na.omit(RootBranch)
  
  D0 = 999*matrix(1,nrow = n,ncol = n)
  D0[RootBranch,] = 0
  D0[,RootBranch] = 0
  
  for (t in 1:length(nonrootlabel)){
    a = nonrootlabel[t]
    Brancha = Branch[,,a]
    Brancha=na.omit(Brancha)
    D0[Brancha,Brancha] = 0
    D0[Brancha,RootBranch] = 0
    D0[RootBranch,Brancha] = 0
  }
  
  D = DISTMAT(X)
  D = D + D0
  
  
  X = t(X)
  
  J <- matrix(nr=nrow(D),nc=ncol(D))
  for (i in 1:nrow(D)){
    J[i,] <- order(D[i,])
  }
  
  
  
  
  skip = 1
  kmin = parameters[1]
  Kmax = parameters[2]
  Kmax = min(Kmax,n)
  eta = parameters[3]
  
  Is<-array(NA, dim = c(1, nrow(label), n))
  
  for (i in 1:n){
    # Contraction
    ki = Kmax
    ratio = rep(0,Kmax)
    while (ki>=kmin){
      Ii = J[i,1:ki]
      Xi = X[,Ii]
      xmean = rowMeans(Xi)
      Xi = t((scale(t(Xi), center = TRUE, scale = FALSE)))
      si = svd(Xi)
      si = si[["d"]]
      ratio[ki] = norm(si[(d+1):length(si)],'2')/norm(si[1:d],'2')
      if (ratio[ki]<eta){
        break
      }else{
        ki = ki-skip
      }
    }
    
    if (ki<kmin){
      j = which.min((ratio[kmin:Kmax]))
      ki =  ki+skip-1+j
      Ii = J[i,1:ki]
      Xi = X[,Ii]
      xmean = rowMeans(Xi)
      Xi = t((scale(t(Xi), center = TRUE, scale = FALSE)))
    }
    Qi = svd(Xi)
    Qi = Qi[["u"]]
    Qi = Qi[,1:d]
    
    # Expansion
    if(ki+1<=Kmax){
      for (j in (ki+1):Kmax){
        xj = X[,J[i,j]]-xmean
        thetaj = t(Qi) %*% xj
        if (norm(xj-Qi %*% thetaj)< eta *norm(thetaj)){
          Ii = c(Ii,J[i,j])
        }
      }
    }
    Is[,1:length(Ii),i]=Ii
    
  }
  
  X = t(X)
  # Compute local information matrix for all datapoints
  print('Compute local information matrices for all datapoints...')
  Bi <- array(NA, dim = c(nrow(label), nrow(label), n))
  for (i in 1:n){
    # Compute correlation matrix W
    Ii = Is[, ,i]
    Ii = na.omit(Ii)
    kt = length(Ii)
    Xi = t(t(X[Ii,])-colMeans(X[Ii,]))
    W = Xi %*% t(Xi)
    W = (W + t(W))/2
    
    # Compute local information by computing d largest eigenvectors of W
    Sch = Schur(W)
    Vi = Sch[["Q"]]
    Si = Sch[["T"]]
    Ji = order(-diag(Si))
    s = (-diag(Si))[order(-diag(Si))]
    if (length(Ji)<d){
      d = length(Ji)
      print(paste("Target dimensionality reduced to",d))
    }
    Vi = as.matrix(Vi[,Ji[1:d]])
    
    if (0){
      # Store eigenvectors in G (Vi is the space with the maximum variance, i.e. a good approximation of the tangent space at point Xi)
      # The constant 1/sqrt(kt) serves as a centering matrix
      Gi = cbind(as.matrix(rep(1/sqrt(kt),kt)),Vi)
      # Compute Bi = I - Gi * Gi'
      bi = diag(rep(1,kt)) - Gi %*% t(Gi)
      Bi[1:nrow(bi),1:ncol(bi) ,i] = bi
    }else{
      Pi = diag(rep(1,kt)) - Vi %*% t(Vi)
      Wi = Pi - matrix(rep(colMeans(Pi),kt*kt),nrow=kt,ncol=kt)
      Bi[1:nrow(Wi),1:ncol(Wi) ,i] = Wi %*% t(Wi)
    }
  }
  
  # Construct sparse matrix B (= alignment matrix)
  print("Construct alignment matrix...")
  B = new("dgTMatrix",
          i = as.integer(1:n-1),
          j = as.integer(1:n-1), x=rep(1,n), Dim=as.integer(c(n,n)))
  for (i in 1:n){
    Ii = Is[, ,i]
    Ii = na.omit(Ii)
    B[Ii,Ii] = B[Ii,Ii] + Bi[1:length(Ii),1:length(Ii),i]   # sum Bi over all points
    B[i,i] = B[i,i] - 1
  }
  B = (B + t(B))/2                                          # make sure B is symmetric
  B = as.matrix(B)
  
  # For sparse datasets, we might end up with NaNs in M. We just set them to zero for now...
  B[B==-Inf]<-0
  B[is.na(B)]<-0
  
  # Perform eigenanalysis of matrix B
  print("Perform eigenanalysis...")
  eigB = eigen(B)
  ind = order(eigB$values, decreasing = FALSE)[1:(d +10)]
  
  # Final embedding coordinates
  D = eigB$values[ind[2:(d + 1)]]
  mappedX = eigB$vectors[, ind[2:(d + 1)]]
  T = mappedX
  
  
  picknonrootbranch = Branch[,,nonrootlabel[1]]
  picknonrootbranch = na.omit(picknonrootbranch)
  nonroot = picknonrootbranch[1]
  
  if (T[root]>T[nonroot]){
    T=-T
  }
  a = min(T)
  b = max(T)
  T = (T-a)/(b-a)
  return(T)
  
}
