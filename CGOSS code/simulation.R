#This R code shows the comparison between different subsample methods on simulated full dataset.
#People can choose 4 types of generated dataset from CASE1-CASE4
#People can also choose 3 different setting of random effect: N.ori, N.ML, N.large, T to generate the dataset.
#People can use "large" or "small" to control the differences between the size of generated groups


#The function Comp is the main function that do the comparison between 5 subsample method. 
#People can change the value of "n" to have different size of subsample set.
#People can change the value of "obj.c" to choose different setting of the objective function.

#We use CPP to help speed up the code. People must put 'myoss.cpp','lmm_fast.cpp' to the environment.




rm(list=ls())
library(devtools)
library(ClusterR)
library(MASS)
Rcpp::sourceCpp('myoss.cpp')
Rcpp::sourceCpp("lmm_fast.cpp")



iboss=function(x,k){
  ind=NULL
  m = ncol(x)
  r = rep(floor(k/2/m),m)
  if(sum(r)<k/2) r[1:((k-2*sum(r))/2)] = r[1:((k-2*sum(r))/2)]+1
  candi=1:nrow(x)
  for (i in 1:m)
  {
    xi = x[candi,i]
    j1 = top_k(xi,r[i])+1
    j2 = bottom_k(xi,r[i])+1
    j = unique(c(j1,j2))
    if(length(j)<2*r[i]) {jj=(1:length(candi))[-j];j=c(j,jj[1:(2*r[i]-length(j))])}
    ind = c(ind,candi[j])
    candi=setdiff(candi,ind)
  }
  return(ind)
}
scalex=function(a){
  2*(a-min(a))/(max(a)-min(a))-1
}
find.sigma=function(xx, yy, beta, nc, R){
  
  eta.hat<-yy-cbind(1, xx)%*%beta
  
  b <- c()
  Ue <- 0
  for (i in 1:length(eta.hat)) {
    b[i] <- sum((eta.hat[i]*rep(1,length(eta.hat)) - eta.hat)^2)/2
    Ue <- Ue + b[i]
  }
  
  if(length(which(nc==0))!=0){
    nc <- nc[-which(nc==0)]
  }
  R.1 <- length(nc)
  
  e1 <- c()
  for(j in 1:nc[1]) {
    e1[j] <- sum((eta.hat[1:nc[1]][j]*rep(1,nc[1]) - eta.hat[1:nc[1]])^2)/(2*nc[1])
  }
  
  
  e <- c()
  e[1] <- sum(e1)
  
  if (R.1 > 1) {  
    for (i in 2:R.1) {
      e2 <- c()
      for (j in 1:nc[i]) {
        e2[j] <- sum((eta.hat[(sum(nc[1:(i-1)]) + j)]*rep(1,nc[i]) - eta.hat[(sum(nc[1:(i-1)]) + 1):sum(nc[1:i])])^2)/(2*nc[i])
      }
      e[i] <- sum(e2)
    }
  }
  
  
  Ua <- sum(e)
  NF <- sum(nc)
  NS <- sum(nc^2)
  Var.e <- Ua/(NF-R)
  if((NF^2-NS)==0){
    Var.a=0
  }else{
    Var.a <- Ue/(NF^2-NS) -  Ua*(NF^2 - NF)/((NF^2 - NS)*(NF-R))
  }
  if (Var.a < 0) {
    Var.a <- 0
  }
  
  
  
  
  return(list(Var.a,Var.e))
}
find.beta=function(xx, yy, Var.a, Var.e, nc, R, p){
  if(length(which(nc==0))!=0){
    nc <- nc[-which(nc==0)]
  }
  R.1 <- length(nc)
  sn <- c(0, cumsum(nc))
  
  
  XVX <- matrix(0, p+1, p+1)
  XVY <- matrix(0, p+1, 1)
  
  for (i in 1:R.1){
    gamma <- (nc[i]*Var.a)/(Var.e+nc[i]*Var.a)
    inv.V <- (1/Var.e)*(diag(1,(nc[i]))-(gamma/(nc[i]))*matrix(1,(nc[i]),(nc[i])))
    if(nc[i] ==1){
      XVX <- XVX + t(cbind(1,t(xx[(sn[i]+1):(sn[i+1]),])))%*%inv.V%*%cbind(1,t(xx[(sn[i]+1):(sn[i+1]),]))
      XVY <- XVY + t(cbind(1,t(xx[(sn[i]+1):(sn[i+1]),])))%*%inv.V%*%yy[(sn[i]+1):(sn[i+1])]
    }else{
      XVX <- XVX + t(cbind(1,xx[(sn[i]+1):(sn[i+1]),]))%*%inv.V%*%cbind(1,xx[(sn[i]+1):(sn[i+1]),])
      XVY <- XVY + t(cbind(1,xx[(sn[i]+1):(sn[i+1]),]))%*%inv.V%*%yy[(sn[i]+1):(sn[i+1])]
    }
  }
  
  
  if (det(XVX) == 0) {
    cat("Nonsingular")
  }
  
  bt<-solve(XVX)%*%XVY
  
  return(bt)
}
Est_hat=function(xx, yy, beta, Var.a, Var.e, nc, R, p){
  
  beta0 <- as.matrix(lm(yy ~ xx)$coefficients)
  
  sigma.hat0 <- find_sigma_cpp(xx, yy, beta0, nc, R)
  
  beta1 <- find_beta_cpp(xx, yy, sigma.hat0[[1]],sigma.hat0[[2]], nc, R, p)
  
  sigma.hat1 <- find_sigma_cpp(xx, yy, beta1, nc, R)
  
  
  beta2 <- find_beta_cpp(xx, yy, sigma.hat1[[1]], sigma.hat1[[2]], nc, R, p)
  
  
  bt.mse <- sum((beta2[-1]-beta)^2)
  sigma.hat2 <- find_sigma_cpp(xx, yy, beta2, nc, R)
  bt0.mse <- (beta2[1] - 1)^2
  va.mse <- sum((sigma.hat2[[1]] - Var.a)^2)
  ve.mse <- sum((sigma.hat2[[2]] - Var.e)^2)
  
  return(list(bt.mse, va.mse, ve.mse, bt0.mse, beta2, sigma.hat2[[1]], sigma.hat2[[2]]))
}
MSPE_fn=function(fy,fx, sx, sy, beta, Var.a, Var.e, nc,C, R){
  index <- 1
  mv_hat <- c()
  for (i in 1:R) {
    mv_hat[i] <- (Var.a/(Var.e+nc[i]*Var.a)) * sum((sy - cbind(1, sx)%*%beta)[index:(index + nc[i] - 1)]) 
    index <- index + nc[i]
  }
  
  
  y_hat <- cbind(1, fx)%*%beta + rep(mv_hat, C)
  mspe <- mean((fy - y_hat)^2)
  
  return(mspe)
}
generate_groups <- function(R, m, N,V) {
  if (N <= R * m) {
    stop("N must be greater than R * m to ensure all integers are greater than m")
  }
  if(V=="large"){Vu<-300*N/(5*R)}
  if(V=="small"){Vu<-N/(5*R)}
  adjusted_sum <- N - R * m
  random_numbers <- abs(rnorm(R, mean = N/R, sd = Vu))
  random_numbers <- random_numbers / sum(random_numbers) * adjusted_sum
  result <- ceiling(random_numbers + m)
  current_sum <- sum(result)
  difference <- N - current_sum
  
  while (difference != 0) {
    if (difference > 0) {
      idx <- sample(1:R, 1)
      result[idx] <- result[idx] + 1
      difference <- difference - 1
    } else {
      idx <- sample(1:R, 1)
      if (result[idx] - 1 > m) {
        result[idx] <- result[idx] - 1
        difference <- difference + 1
      }
    }
  }
  
  return(result)
}
assign_clusters <- function(data, centroids) {
  cluster_assignments <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    distances <- apply(centroids, 1, function(centroid) sum((data[i, ] - centroid) ^ 2))
    cluster_assignments[i] <- which.min(distances)
  }
  
  return(cluster_assignments)
}
mbky <- function(setseed, FXX, y, n, Cn) {
  set.seed(setseed)
  
  repeat {
    sz=550 - 100 * log10(nrow(FXX) + 1)
    mini_batch_kmeans <- MiniBatchKmeans(FXX, clusters = Cn, batch_size = sz, num_init = 5, max_iters = 20, initializer = 'kmeans++')
    centroids <- mini_batch_kmeans$centroids
    batchs <- assign_clusters(FXX, centroids)
    cluster_sizes <- table(batchs)
    
    threshold <- n / Cn
    if (any(cluster_sizes < threshold)) {
      Cn <- Cn - 1  
    } else {
      break  
    }
  }
  
  R_CGOSS = length(cluster_sizes)
  original_indices <- 1:nrow(FXX)
  data_with_cluster <- data.frame(FXX, y = y, cluster = batchs, original_index = original_indices)
  data_sorted <- data_with_cluster[order(data_with_cluster$cluster), ]
  data_matrix_sorted <- as.matrix(data_sorted[, !(names(data_sorted) %in% c("row.names", "cluster", "original_index", "y")), drop = FALSE])
  sorted_y <- data_sorted$y
  cluster_sizes <- table(data_sorted$cluster)
  cluster_sizes_vector <- as.vector(cluster_sizes)
  sorted_indices <- data_sorted$original_index
  return(list(R_CGOSS = R_CGOSS, 
              data_matrix_sorted = data_matrix_sorted, 
              sorted_y = sorted_y, 
              cluster_sizes_vector = cluster_sizes_vector, 
              sorted_indices = sorted_indices))
}
count.info=function(xx,yy,nc,R,p){
  xx <- apply(xx, 2, scalex)
  
  if(length(nc)==1){
    X <- model.matrix(~ xx)
    model <- lm(yy ~ xx)
    residuals <- resid(model)
    RSS <- sum(residuals^2)
    Var.lm<-RSS/(nrow(X)-p-1)
    D<-det(t(X) %*% X)/(Var.lm^p)
    A <- Var.lm*sum(diag(solve( t(X) %*% X)) )
  }
  else{
    
    
    
    beta0 <- as.matrix(lm(yy ~ xx)$coefficients)
    sigma.hat0 <- find_sigma_cpp(xx, yy, beta0, nc, R)
    beta1 <- find_beta_cpp(xx, yy, sigma.hat0[[1]],sigma.hat0[[2]], nc, R, p)
    sigma.hat1 <- find_sigma_cpp(xx, yy, beta1, nc, R)
    beta2 <- find_beta_cpp(xx, yy, sigma.hat1[[1]], sigma.hat1[[2]], nc, R, p)
    sigma.hat2 <- find_sigma_cpp(xx, yy, beta2, nc, R)
    
    Var.a<-sigma.hat2[[1]]
    Var.e<-sigma.hat2[[2]]
    if(length(which(nc==0))!=0){
      nc <- nc[-which(nc==0)]
    }
    R.1 <- length(nc)
    sn <- c(0, cumsum(nc))
    
    
    XVX <- matrix(0, p+1, p+1)
    for (i in 1:R.1){
      
      gamma <- (nc[i]*Var.a)/(Var.e+nc[i]*Var.a)
      inv.V <- (1/Var.e)*(diag(1,(nc[i]))-(gamma/(nc[i]))*matrix(1,(nc[i]),(nc[i])))
      
      if (any(is.nan(inv.V)) || any(is.infinite(inv.V))) {
        cat("inv.V  NaN or
            Inf\n")
      }
      
      if(nc[i] ==1){
        XVX <- XVX + t(cbind(1,t(xx[(sn[i]+1):(sn[i+1]),])))%*%inv.V%*%cbind(1,t(xx[(sn[i]+1):(sn[i+1]),]))
      }else{
        XVX <- XVX + t(cbind(1,xx[(sn[i]+1):(sn[i+1]),]))%*%inv.V%*%cbind(1,xx[(sn[i]+1):(sn[i+1]),])
      }
      
    }
    
    D<-det(XVX)
    A<-sum(diag( solve(XVX) ))
  }
  
  
  return(c(D,A))
}
findsubforCGOSS<-function(n,R){
  if (n %% R != 0) {
    me=floor(n/R)
    loss=n-me*R
    
    mCGOSS=c(rep(me+1,loss),rep(me,R-loss))
  }
  else{
    mCGOSS=c(rep(n/R,R))
  }
  return(mCGOSS)
}
GOSS<-function(setseed,FXX,FY,n,Cn,p){
  cluster=mbky(setseed,FXX,FY,n,Cn)
  R_CGOSS=cluster$R_CGOSS
  FXXXX <- cluster$data_matrix_sorted
  FYYY<- cluster$sorted_y
  SCC <- c(0, cumsum(cluster$cluster_sizes_vector))
  mcgoss<-findsubforCGOSS(n,R_CGOSS)
  index.CGOSS <- integer(0) 
  for (i in 1:(length(SCC) - 1)) {
    current_indices <- (SCC[i] + 1):SCC[i + 1]
    index.CGOSS <- c(index.CGOSS, OAJ2_cpp(apply(FXXXX[current_indices, ], 2, scalex), mcgoss[i], tPow=2) + SCC[i])
  }
  index_CGOSS_interation <- cluster$sorted_indices[index.CGOSS]
  ncCGOSS <- mcgoss
  D.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)[1]
  A.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)[2]
  return(list(index = index_CGOSS_interation,D = D.after,A = A.after,R = R_CGOSS,nc = ncCGOSS,C=cluster$cluster_sizes_vector))
}

Comp=function(N_all,p, R, Var.e, nloop, n, dist_x="case1", dist_a="NaN",groupsize,anneal){
  big_column_vector<-c()
  beta=rep(1, p)
  m=ceiling(n / R)
  sigma=diag(0.5,p,p)+matrix(0.5,p,p)
  lrs=length(N_all)
  names=c("CGOSS.bt.mat",  "IBOSS.bt.mat", "OSS.bt.mat", 
          "knowGOSS.bt.mat", "GOSS.bt.mat","GIBOSS.bt.mat","knowGIBOSS.bt.mat",
          "CGOSS.pred",  "IBOSS.pred", "OSS.pred", 
          "knowGOSS.pred", "GOSS.pred","GIBOSS.pred","knowGIBOSS.pred",
          "CGOSS.bt0.dif","IBOSS.bt0.dif","OSS.bt0.dif",
          "knowGOSS.bt0.dif","GOSS.bt0.dif","GIBOSS.bt0.dif","knowGIBOSS.bt0.dif",
          "CGOSS.Var.a","IBOSS.Var.a","OSS.Var.a",
          "knowGOSS.Var.a","GOSS.Var.a","GIBOSS.Var.a","knowGIBOSS.Var.a",
          "CGOSS.Var.e","IBOSS.Var.e","OSS.Var.e",
          "knowGOSS.Var.e","GOSS.Var.e","GIBOSS.Var.e","knowGIBOSS.Var.e")
  mat_names=c("CGOSS.bt",  "IBOSS.bt", "OSS.bt", 
              "knowGOSS.bt", "GOSS.bt","GIBOSS.bt","knowGIBOSS.bt")
  for(name in names) {
    assign(name, matrix(NA, 1, nloop*lrs), envir = .GlobalEnv)
  }
  for(name in mat_names) {
    assign(name, matrix(NA, p+1, nloop*lrs), envir = .GlobalEnv)
  }
  itr = 0
  #######
  for (j in 1:lrs) {
    D.CGOSS=0
    A.CGOSS=0
    A.OSS=0
    time.CGOSS=0
    meanR=0
    for (k in 1:nloop) {
      
      N<-N_all[j]
      random_numbers <- generate_groups(R,m,N,groupsize)
      C <- round(random_numbers)
      SC = c(0, cumsum(C))
      if (k%/%100 == k/100) cat(k, "-")
      itr <- itr+1
      set.seed(k* 100000)
      if(dist_a == "N.ori") {Var.a = 0.5; Fa = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C)}
      if(dist_a == "N.ML") { Var.a <- 0; Fa <- 0 }
      if(dist_a == "N.large") {Var.a = 100; Fa = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C)}
      if(dist_a=="T"){Var.a = 3;Fa = rep(rt(R,3), C)}
      
      Fe = rnorm(max(SC),mean = 0,sd = sqrt(Var.e))
      FXX = matrix(0, nrow = max(SC), ncol = p)
      index.knowGOSS <-index.CGOSS<- index.GOSS<- index.GIBOSS<- index.knowGIBOSS <- c()
      cpu_time_index_goss<-0
      
      
      ##############
      for (i in 1:R) {
        
        setseed =  k * 100000 + i * 100
        set.seed(setseed)
        if(dist_x=="case1") {FXX[(SC[i] + 1):(SC[i+1]),]=matrix(runif(C[i]*p, -1, 1),C[i],p)}
        if(dist_x=="case2") {FXX[(SC[i] + 1):(SC[i+1]),]=mvrnorm(C[i], rep(0, p), sigma)}
        if(dist_x=="case3") {FXX[(SC[i] + 1):(SC[i+1]),]=matrix(runif(C[i]*p, -1.55+i/20, 0.45+i/20),C[i],p)}
        if(dist_x=="case4") {FXX[(SC[i] + 1):(SC[i+1]),]=mvrnorm(C[i], rep(-2+(i-1)/5, p), sigma) }
        index.knowGOSS <- c(index.knowGOSS, OAJ2_cpp(apply(FXX[(SC[i] + 1):(SC[i+1]),],2,scalex),m, tPow=2) + SC[i])
        index.knowGIBOSS <- c(index.knowGIBOSS, iboss(FXX[(SC[i] + 1):(SC[i+1]),],m) + SC[i])
        Fori <- FXX
      }
      FYori <- 1 + Fori%*%beta + Fa + Fe
      shuffled_indices <- sample(nrow(FXX))
      shuffled_df <- FXX[shuffled_indices, ]
      rownames(shuffled_df) <- rownames(FXX)[shuffled_indices]
      FXX <- shuffled_df
      CC=rep(N/R,R)
      SSC = c(0, cumsum(CC))
      for(i in 1:R){
        
        index.GOSS <- c(index.GOSS, OAJ2_cpp(apply(FXX[(SSC[i] + 1):(SSC[i+1]),],2,scalex),m, tPow=2) + SSC[i])
        
        index.GIBOSS <- c(index.GIBOSS, iboss(FXX[(SSC[i] + 1):(SSC[i+1]),],m) + SSC[i])
      }

      FY <- 1 + FXX%*%beta + Fa + Fe
      nc <- rep(m,R)
      
      ########################################## OSS
      nc2 <- c()
      index.OSS <- sort(OAJ2_cpp(apply(FXX,2,scalex),n, tPow=2))
      for (i in 1:R) {
        nc.OSS <- which(index.OSS >= (SSC[i] + 1) & index.OSS <= (SSC[i+1]))
        nc2[i] <- length(nc.OSS)
      }
      
      #########################CGOSS###########################################################################
      
      time.start<-Sys.time()
      D.oss=count_info_cpp(FXX[index.OSS,],FY[index.OSS,],n,1,p)[1]
      A.oss=count_info_cpp(FXX[index.OSS,],FY[index.OSS,],n,1,p)[2]
      obj.best<-obj.all<-(anneal/p)*log(D.oss)-(1-anneal)*(log(A.oss))
      obj.before<-obj.best
      index.best<-index.OSS
      nc.best<-n
      R.best<-1
      Cn=2
      C.CGOSS=C
      T.initial<-n*p/5
      repeat {
        informat <- GOSS(setseed, FXX, FY, n, Cn, p)
        obj.candi <- (anneal/p)*log(informat$D) - (1-anneal)*(log(informat$A/p))
        if (obj.candi >= obj.best) {
          if (obj.candi == obj.best && obj.best == obj.before) {
            break
          }
          obj.best    <- obj.candi
          index.best  <- informat$index
          nc.best     <- informat$nc
          C.CGOSS     <- informat$C
          R.best      <- informat$R
          obj.before  <- obj.candi
          Cn <- Cn + 1
          next
        }
        alpha <- if (obj.before == obj.candi) 0.3 else if (Cn == informat$R) 0.95 else 0.85
        T.cool <- T.initial * alpha^Cn
        heatprob <- exp(-(obj.best - obj.candi)/T.cool)
        if (runif(1) < heatprob) {
          Cn <- Cn + 1
          obj.before <- obj.candi
          next
        } else {
          break
        }
      }
      
      final_index_CGOSS<-index.best
      ncCGOSS<-nc.best
      R_CGOSS<-R.best
      time.end<-Sys.time()
      time.CGOSS<-time.CGOSS+as.numeric(difftime(time.end, time.start, units = "secs"))
      print(time.CGOSS)
      print(R_CGOSS)
      meanR<-meanR+R_CGOSS

      ############################################## IBOSS
      nc3 <- c()
      index.IBOSS <- sort(iboss(FXX,n))
      for (i in 1:R) {
        nc.IBOSS <- which(index.IBOSS >= (SC[i] + 1) & index.IBOSS <= (SC[i+1]))
        nc3[i] <- length(nc.IBOSS)
      }
      
      
      ##########################################################  GOSS

      GOSS.Est <- Est_hat_cpp(xx=FXX[index.GOSS,], yy=FY[index.GOSS,], 
                              beta, Var.a, Var.e, nc, R, p)
      GOSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.GOSS,], FY[index.GOSS,], 
                                 GOSS.Est[[5]], GOSS.Est[[6]], GOSS.Est[[7]], nc,C, R)
      GOSS.bt.mat[,itr] <- GOSS.Est[[1]]
      GOSS.Var.a[,itr]<- GOSS.Est[[2]]
      GOSS.Var.e[,itr]<- GOSS.Est[[3]]
      GOSS.bt0.dif[,itr] <- GOSS.Est[[4]]
      GOSS.bt[,itr] <- GOSS.Est[[5]]
      
      ############################################################# estimate knowGOSS
      
      knowGOSS.Est <- Est_hat_cpp(xx=Fori[index.knowGOSS,], yy=FYori[index.knowGOSS,], 
                                  beta, Var.a, Var.e, nc, R, p)
      knowGOSS.pred[,itr] <- MSPE_fn(FYori, Fori, Fori[index.knowGOSS,], FYori[index.knowGOSS,], 
                                     knowGOSS.Est[[5]], knowGOSS.Est[[6]], knowGOSS.Est[[7]], nc,C, R)
      knowGOSS.bt.mat[,itr] <- knowGOSS.Est[[1]]
      knowGOSS.Var.a[,itr]<- knowGOSS.Est[[2]]
      knowGOSS.Var.e[,itr]<- knowGOSS.Est[[3]]
      knowGOSS.bt0.dif[,itr] <- knowGOSS.Est[[4]]
      knowGOSS.bt[,itr] <- knowGOSS.Est[[5]]
      
      ############################################################# estimate GIBOSS
      GIBOSS.Est <- Est_hat_cpp(xx=FXX[index.GIBOSS,], yy=FY[index.GIBOSS,], 
                                beta, Var.a, Var.e, nc, R, p)
      GIBOSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.GIBOSS,], FY[index.GIBOSS,], 
                                   GIBOSS.Est[[5]], GIBOSS.Est[[6]], GIBOSS.Est[[7]], nc,C, R)
      GIBOSS.bt.mat[,itr] <- GIBOSS.Est[[1]]
      GIBOSS.Var.a[,itr]<- GIBOSS.Est[[2]]
      GIBOSS.Var.e[,itr]<- GIBOSS.Est[[3]]
      GIBOSS.bt0.dif[,itr] <- GIBOSS.Est[[4]]
      GIBOSS.bt[,itr] <- GIBOSS.Est[[5]]
      
      ############################################################# estimate knowGIBOSS
      knowGIBOSS.Est <- Est_hat_cpp(xx=Fori[index.knowGIBOSS,], yy=FYori[index.knowGIBOSS,], 
                                    beta, Var.a, Var.e, nc, R, p)
      knowGIBOSS.pred[,itr] <- MSPE_fn(FYori, Fori, Fori[index.knowGIBOSS,], FYori[index.knowGIBOSS,], 
                                       knowGIBOSS.Est[[5]], knowGIBOSS.Est[[6]], knowGIBOSS.Est[[7]], nc,C, R)
      knowGIBOSS.bt.mat[,itr] <- knowGIBOSS.Est[[1]]
      knowGIBOSS.Var.a[,itr]<- knowGIBOSS.Est[[2]]
      knowGIBOSS.Var.e[,itr]<- knowGIBOSS.Est[[3]]
      knowGIBOSS.bt0.dif[,itr] <- knowGIBOSS.Est[[4]]
      knowGIBOSS.bt[,itr] <- knowGIBOSS.Est[[5]]
      ############################################################# estimate CGOSS
      CGOSS.Est <- Est_hat_cpp(xx=FXX[final_index_CGOSS,], yy=FY[final_index_CGOSS,], 
                               beta, Var.a, Var.e, ncCGOSS, R_CGOSS, p)
      CGOSS.pred[,itr]  <- MSPE_fn(FY,FXX , FXX[final_index_CGOSS,], FY[final_index_CGOSS,], 
                                   CGOSS.Est[[5]], CGOSS.Est[[6]], CGOSS.Est[[7]], ncCGOSS,C.CGOSS, R_CGOSS)
      CGOSS.bt.mat[,itr] <- CGOSS.Est[[1]]
      CGOSS.Var.a[,itr]<- CGOSS.Est[[2]]
      CGOSS.Var.e[,itr]<- CGOSS.Est[[3]]
      CGOSS.bt0.dif[,itr] <- CGOSS.Est[[4]]
      CGOSS.bt[,itr] <- CGOSS.Est[[5]]
      ##########################################################  OSS
      OSS.Est <- Est_hat_cpp(xx=FXX[index.OSS,], yy=FY[index.OSS,], 
                             beta, Var.a, Var.e, nc2, R, p)
      OSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.OSS,], FY[index.OSS,], 
                                OSS.Est[[5]], OSS.Est[[6]], OSS.Est[[7]], nc2,(N/n)*nc2, R)
      OSS.bt.mat[,itr] <- OSS.Est[[1]]
      OSS.Var.a[,itr]<- OSS.Est[[2]]
      OSS.Var.e[,itr]<- OSS.Est[[3]]
      OSS.bt0.dif[,itr] <- OSS.Est[[4]]
      OSS.bt[,itr] <- OSS.Est[[5]]
      ############################################################# estimate IBOSS
      IBOSS.Est <- Est_hat_cpp(xx=FXX[index.IBOSS,], yy=FY[index.IBOSS,], 
                               beta, Var.a, Var.e, nc3, R, p)
      IBOSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.IBOSS,], FY[index.IBOSS,], 
                                  IBOSS.Est[[5]], IBOSS.Est[[6]], IBOSS.Est[[7]], nc3,(N/n)*nc3, R)
      IBOSS.bt.mat[,itr] <- IBOSS.Est[[1]]
      IBOSS.Var.a[,itr]<- IBOSS.Est[[2]]
      IBOSS.Var.e[,itr]<- IBOSS.Est[[3]]
      IBOSS.bt0.dif[,itr] <- IBOSS.Est[[4]]
      IBOSS.bt[,itr] <- IBOSS.Est[[5]]
      cat(j,"-",k,"\n")
    }
    
    
    cat("mean time is",time.CGOSS/nloop,"\n")
    cat("mean CGOSS R is",meanR/nloop,"\n")
    
    
    cat("\n\n")
  }
  
  
  ##########################################################
  mse.Goss<- mse.knowGOSS <- mse.CGOSS <- mse.oss <- mse.iboss <-mse.GIBOSS<-mse.knowGIBOSS<- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.knowGOSS <- c(mse.knowGOSS, mean(knowGOSS.bt.mat[,loc]))
    mse.Goss <- c(mse.Goss, mean(GOSS.bt.mat[,loc]))
    mse.CGOSS <- c(mse.CGOSS, mean(CGOSS.bt.mat[,loc]))
    mse.iboss <- c(mse.iboss, mean(IBOSS.bt.mat[,loc]))
    mse.oss <- c(mse.oss, mean(OSS.bt.mat[,loc]))
    mse.GIBOSS <- c(mse.GIBOSS, mean(GIBOSS.bt.mat[,loc]))
    mse.knowGIBOSS <- c(mse.knowGIBOSS, mean(knowGIBOSS.bt.mat[,loc]))
    
  }
  rec1<-cbind(mse.CGOSS, mse.iboss, mse.oss, mse.knowGOSS, mse.Goss,mse.GIBOSS,mse.knowGIBOSS)
  ###############################################
  mse.CGOSS.Var.a <- mse.oss.Var.a <- mse.iboss.Var.a  <- mse.Goss.Var.a <- mse.GIBOSS.Var.a<- mse.knowGIBOSS.Var.a<- mse.knowGOSS.Var.a <-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    mse.CGOSS.Var.a <- c(mse.CGOSS.Var.a, mean(CGOSS.Var.a[,loc]))
    mse.iboss.Var.a <- c(mse.iboss.Var.a, mean(IBOSS.Var.a[,loc]))
    mse.oss.Var.a <- c(mse.oss.Var.a, mean(OSS.Var.a[,loc]))
    mse.knowGOSS.Var.a <- c(mse.knowGOSS.Var.a, mean(knowGOSS.Var.a[,loc]))
    mse.Goss.Var.a <- c(mse.Goss.Var.a, mean(GOSS.Var.a[,loc]))
    mse.GIBOSS.Var.a <- c(mse.GIBOSS.Var.a, mean(GIBOSS.Var.a[,loc]))
    mse.knowGIBOSS.Var.a <- c(mse.knowGIBOSS.Var.a, mean(knowGIBOSS.Var.a[,loc]))
    
  }
  rec2<-cbind(mse.CGOSS.Var.a, mse.iboss.Var.a, mse.oss.Var.a, mse.knowGOSS.Var.a, mse.Goss.Var.a,mse.GIBOSS.Var.a,mse.knowGIBOSS.Var.a)
  ###############################################
  mse.CGOSS.Var.e <- mse.oss.Var.e <- mse.iboss.Var.e <- mse.Goss.Var.e<- mse.GIBOSS.Var.e<- mse.knowGIBOSS.Var.e <- mse.knowGOSS.Var.e <-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    mse.CGOSS.Var.e <- c(mse.CGOSS.Var.e, mean(CGOSS.Var.e[,loc]))
    mse.iboss.Var.e <- c(mse.iboss.Var.e, mean(IBOSS.Var.e[,loc]))
    mse.oss.Var.e <- c(mse.oss.Var.e, mean(OSS.Var.e[,loc]))
    mse.knowGOSS.Var.e <- c(mse.knowGOSS.Var.e, mean(knowGOSS.Var.e[,loc]))
    mse.Goss.Var.e <- c(mse.Goss.Var.e, mean(GOSS.Var.e[,loc]))
    mse.GIBOSS.Var.e <- c(mse.GIBOSS.Var.e, mean(GIBOSS.Var.e[,loc]))
    mse.knowGIBOSS.Var.e <- c(mse.knowGIBOSS.Var.e, mean(knowGIBOSS.Var.e[,loc]))
    
  }
  rec3<-cbind(mse.CGOSS.Var.e, mse.iboss.Var.e, mse.oss.Var.e, mse.knowGOSS.Var.e, mse.Goss.Var.e,mse.GIBOSS.Var.e,mse.knowGIBOSS.Var.e)
  ##################################################
  mspe.Goss <- mspe.knowGOSS <- mspe.CGOSS<- mspe.oss<- mspe.GIBOSS<- mspe.knowGIBOSS<- mspe.iboss <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    mspe.knowGOSS <- c(mspe.knowGOSS, mean(knowGOSS.pred[,loc]))
    mspe.Goss <- c(mspe.Goss, mean(GOSS.pred[,loc]))
    mspe.CGOSS <- c(mspe.CGOSS, mean(CGOSS.pred[,loc]))
    mspe.iboss <- c(mspe.iboss, mean(IBOSS.pred[,loc]))
    mspe.oss <- c(mspe.oss, mean(OSS.pred[,loc]))
    mspe.GIBOSS <- c(mspe.GIBOSS, mean(GIBOSS.pred[,loc]))
    mspe.knowGIBOSS <- c(mspe.knowGIBOSS, mean(knowGIBOSS.pred[,loc]))
    
  }
  
  rec4 <- cbind(mspe.CGOSS, mspe.iboss, mspe.oss, mspe.knowGOSS, mspe.Goss,mspe.GIBOSS,mspe.knowGIBOSS)
  
  ################################################
  mse.bt0.Goss <- mse.bt0.knowGOSS <- mse.bt0.CGOSS<- mse.bt0.oss<- mse.bt0.GIBOSS<- mse.bt0.knowGIBOSS<- mse.bt0.iboss <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    
    mse.bt0.knowGOSS <- c(mse.bt0.knowGOSS, mean(knowGOSS.bt0.dif[,loc]))
    mse.bt0.Goss <- c(mse.bt0.Goss, mean(GOSS.bt0.dif[,loc]))
    mse.bt0.CGOSS <- c(mse.bt0.CGOSS, mean(CGOSS.bt0.dif[,loc]))
    mse.bt0.iboss <- c(mse.bt0.iboss, mean(IBOSS.bt0.dif[,loc]))
    mse.bt0.oss <- c(mse.bt0.oss, mean(OSS.bt0.dif[,loc]))
    mse.bt0.GIBOSS <- c(mse.bt0.GIBOSS, mean(GIBOSS.bt0.dif[,loc]))
    mse.bt0.knowGIBOSS <- c(mse.bt0.knowGIBOSS, mean(knowGIBOSS.bt0.dif[,loc]))
    
  }
  
  rec5 <- cbind(mse.bt0.CGOSS, mse.bt0.iboss, mse.bt0.oss, mse.bt0.knowGOSS, mse.bt0.Goss,mse.bt0.GIBOSS,mse.bt0.knowGIBOSS)
  
  save(rec1, rec2, rec3, rec4, rec5, file = paste0(dist_a,"_", dist_x,".Rdata"))
  
  return(list(rec1,rec2,rec3,rec4,rec5))
}

#########################



N=c(1e4,5e4,1e5)

result = Comp(N,p=50,R=20,Var.e=9,nloop=200,n=1e3,dist_x ="case4", dist_a="N.ori",groupsize="large",anneal=0.1)
result


plot(log10(N), log10(result[[1]][,1]), xlab=expression(log["10"](N)), ylab = expression(log["10"](MSE)), pch=1,lwd=2,
     ylim=c(min(log10(result[[1]])), max(log10(result[[1]]))),xlim=c(min(log10(N)), max(log10(N))), type="o",main = "")
pp <- dim(result[[1]])[2]
for(i in 2:pp){
  lines(log10(N), log10(result[[1]][,i]), type="o", pch=(i), lty=i, col=i,lwd=2)
}

legend("topright", lty=1:pp, pch=(1:pp), col=1:pp,
       legend=c("CGOSS","IBOSS","OSS","known-GOSS","GOSS","GIBOSS","knowGIBOSS"),cex=0.75)


###########################################

plot(log10(N), log10(result[[4]][,1]), xlab=expression(log["10"](N)), ylab = expression(log["10"](MSPE)), pch=1,lwd=2,
     ylim=c(min(log10(result[[4]])), max(log10(result[[4]]))),xlim=c(min(log10(N)), max(log10(N))), type="o",main = "")
pp <- dim(result[[4]])[2]
for(i in 2:pp){
  lines(log10(N), log10(result[[4]][,i]), type="o", pch=(i), lty=i, col=i,lwd=2)
}

legend("topright", lty=1:pp, pch=(1:pp), col=1:pp,
       legend=c("CGOSS","IBOSS","OSS","known-GOSS","GOSS","GIBOSS","knowGIBOSS"),cex=0.75)




















