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
  
  Is<-array(list(), n)
  
  
  for (i in 1:n){
    # Contraction
    ki = Kmax
    ratio = rep(0,Kmax)
    while (ki>=kmin){
      Ii = J[i,1:ki]
      Xi = X[,Ii]
      xmean = rowMeans(Xi)
      Xi =  Xi - matrix(rep(xmean, ki), ncol = ki)
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
      Xi =  Xi - matrix(rep(xmean, ki), ncol = ki)
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
    Is[[i]]=Ii
    
  }
  
  
  X = t(X)
  # Compute local information matrix for all datapoints
  print('Compute local information matrices for all datapoints...')
  Bi = array(list(), n)
  K = numeric(n)
  for (i in 1:n) {
    K[i] = length(Is[[i]])
  }
  
  for (i in 1:n){
    # Compute correlation matrix W
    Ii = Is[[i]]
    kt = K[i]
    Xi = X[Ii,] - matrix(rep(colMeans(X[Ii, ]), kt), nrow = kt, byrow = TRUE)
    W = Xi %*% t(Xi)
    W = (W + t(W))/2
    
    # Compute local information by computing d largest eigenvectors of W
    tmp = eigs_sym(W, d, which = "LM")
    Vi = tmp$vectors[,1:d]
    
    Pi = diag(rep(1,kt)) - Vi %*% t(Vi)
    Wi = Pi - matrix(rep(colMeans(Pi), kt), nrow = kt, byrow = TRUE)
    Bi[[i]] = Wi %*% t(Wi)
    
  }
  
  
  # Construct sparse matrix B (= alignment matrix)
  print("Construct alignment matrix...")
  B = diag(rep(1, n))
  for (i in 1:n){
    Ii = Is[[i]]
    B[Ii,Ii] = B[Ii,Ii] + Bi[[i]]  # sum Bi over all points
    B[i,i] = B[i,i] - 1
    
  }
  B = (B + t(B))/2                                          # make sure B is symmetric
  
  
  # For sparse datasets, we might end up with NaNs in M. We just set them to zero for now...
  B[B==-Inf]<-0
  B[is.na(B)]<-0
  
  B <- as(B,'dgCMatrix')
  
  # Perform eigenanalysis of matrix B
  print("Perform eigenanalysis...")
  
  eigB = eigs_sym(B, d+10, which = "SM")
  ind = order(eigB$values, decreasing = FALSE)
  
  # Final embedding coordinates
  D = eigB$values[ind[2:(d + 1)]]
  mappedX = eigB$vectors[, ind[2:(d + 1)]]
  T = mappedX
  
  
  
  picknonrootbranch = Branch[,,nonrootlabel[1]]
  picknonrootbranch = na.omit(picknonrootbranch)
  nonroot = picknonrootbranch[1]
  
  a = min(T)
  b = max(T)
  T = (T-a)/(b-a)
  
  id = which(T<T[root])
  T[id] = 2*T[root]-T[id]
  
  
  pseudo = T
  return(pseudo)
  
}
