#This R code shows the comparison between different subsample methods on a Realdata Set: "climate_change_impact_on_agriculture_2024.csv"

#We choose "day" as the groups in this dataset but people can try other settings by change the grouping option. 
#For different setting we recommend to set different step length for CGOSS to ensure the calculation speed quick enough.

#The function Comp is the main function that do the comparison between 5 subsample method. 
#People can change the value of "n" to have different size of subsample set.
#People can change the value of "obj.c" to choose different setting of the objective function.

#We use CPP to help speed up the code. People must put 'myoss.cpp','lmm_fast.cpp' to the environment.



#######################
rm(list=ls())

library(devtools)
library(MASS)
library(ClusterR)

Rcpp::sourceCpp('myoss.cpp')
Rcpp::sourceCpp('lmm_fast.cpp')

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
MSPE_fn=function(fy,fx, sx, sy, beta, Var.a, Var.e, nc,C, R){
  index <- 1
  mv_hat <- c()
  for (i in 1:R) {
    mv_hat[i] <- (Var.a/(Var.e+nc[i]*Var.a)) * sum((sy - cbind(1, sx)%*%beta)[index:(index + nc[i] - 1)]) 
    index <- index + nc[i]
  }
  
  y_hat <- cbind(1, fx)%*%beta - rep(mv_hat, C)
  mspe <- mean((fy - y_hat)^2)
  
  return(mspe)
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
   
    mini_batch_kmeans <- MiniBatchKmeans(FXX, clusters = Cn, batch_size = n, num_init = 3, max_iters = 5, initializer = 'kmeans++')
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
  return(list(index = index_CGOSS_interation,D = D.after,A = A.after,R = R_CGOSS,C=cluster$cluster_sizes_vector
              ,FX=FXXXX,FY=FYYY,nc = ncCGOSS))
}
MSE_LM<-function(xx,yy,beta){
  p<-ncol(xx)
  beta0 <- as.matrix(lm(yy ~ xx)$coefficients)
  mse<-sum((beta0[-1]-beta)^2)
  bt0<-(beta0[1]-1)^2
  return(list(bt0,mse,beta0))
}
MSPE_LM<-function(xx,yy,beta){
  n<-nrow(xx)
  y.est<-cbind(1,xx)%*%beta
  mspe<- mean((yy-y.est)^2)
}
Comp=function(X,Y,SSC,groupsize,n,setted_cluster){
  R=length(SSC)-1
  p=ncol(X)
  N=nrow(X)
  
  names=c("CGOSS.bt.mat",  "IBOSS.bt.mat", "OSS.bt.mat", 
          "GOSS.bt.mat","GIBOSS.bt.mat",
          "CGOSS.pred",  "IBOSS.pred", "OSS.pred", 
          "GOSS.pred","GIBOSS.pred",
          "CGOSS.bt0.dif","IBOSS.bt0.dif","OSS.bt0.dif",
          "GOSS.bt0.dif","GIBOSS.bt0.dif",
          "CGOSS.Var.a","IBOSS.Var.a","OSS.Var.a",
          "GOSS.Var.a","GIBOSS.Var.a",
          "CGOSS.Var.e","IBOSS.Var.e","OSS.Var.e",
          "GOSS.Var.e","GIBOSS.Var.e")
  mat_names=c( "CGOSS.bt",  "IBOSS.bt", "OSS.bt", 
               "GOSS.bt","GIBOSS.bt")
  for(name in names) {
    assign(name, matrix(NA, 1, 1), envir = .GlobalEnv)
  }
  for(name in mat_names) {
    assign(name, matrix(NA, p+1, 1), envir = .GlobalEnv)
  }
  FXX<-X
  FY <- matrix(Y, ncol = 1)
  C<-groupsize
  
  
  time.start<-Sys.time()
  beta0 <- as.matrix(lm(FY ~ FXX)$coefficients)
  sigma.hat0 <- find_sigma_cpp(FXX, FY, beta0, C, R)
  beta1 <- find_beta_cpp(FXX, FY, sigma.hat0[[1]],sigma.hat0[[2]], C, R, p)
  sigma.hat1 <- find_sigma_cpp(FXX, FY, beta1, C, R)
  beta.2 <- find_beta_cpp(FXX, FY, sigma.hat1[[1]], sigma.hat1[[2]], C, R, p)
  beta=beta.2[-1]
  sigma.hat2 <- find_sigma_cpp(FXX, FY, beta.2, C, R)
  Var.a=sigma.hat2[[1]]
  Var.e=sigma.hat2[[2]] 
  
  time.end<-Sys.time()
  time.full<-as.numeric(difftime(time.end, time.start, units = "secs"))
  
  print("get full_beta")
  time.CGOSS=0
  setseed=set.seed(42)
  index.GOSS<- index.GIBOSS<- c()
  
  
  m <- ceiling(n / R)
  nc<- findsubforCGOSS(n,R)
  nc.GIBOSS<-c()
  for(i in 1:R){
    index.GOSS <- c(index.GOSS, OAJ2_cpp(apply(FXX[(SSC[i] + 1):(SSC[i+1]),],2,scalex),nc[i], tPow=2) + SSC[i])
    index.GIBOSS <- c(index.GIBOSS, iboss(FXX[(SSC[i] + 1):(SSC[i+1]),],nc[i]) + SSC[i])
    nc.GIBOSS<-c(nc.GIBOSS,length(iboss(FXX[(SSC[i] + 1):(SSC[i+1]),],nc[i])))
  }
  ########################################## OSS
  nc2 <- c()
  index.OSS <- sort(OAJ2_cpp(apply(FXX,2,scalex),n, tPow=2))
  for (i in 1:R) {nc.OSS <- which(index.OSS >= (SSC[i] + 1) & index.OSS <= (SSC[i+1]))
  nc2[i] <- length(nc.OSS)
  }
  ############################################## IBOSS
  nc3 <- c()
  index.IBOSS <- sort(iboss(FXX,n))
  for (i in 1:R) {nc.IBOSS <- which(index.IBOSS >= (SSC[i] + 1) & index.IBOSS <= (SSC[i+1]))
  nc3[i] <- length(nc.IBOSS)
  }
  ############################################## CGOSS   
  
  time2.start<-Sys.time()
  
  informat <- GOSS(setseed, FXX, FY, n, setted_cluster, p)
  final_index_CGOSS  <- informat$index
  ncCGOSS     <- informat$nc
  C.est       <- informat$C
  FX.est      <- informat$FX
  FY.est      <- informat$FY
  R_CGOSS      <- informat$R
  
  time2.end<-Sys.time()
  time.CGOSS<-time.CGOSS+as.numeric(difftime(time2.end, time2.start, units = "secs"))
  print(time.CGOSS)
  print(R_CGOSS)
  ##########################################################  GOSS
  
  GOSS.Est <- Est_hat_cpp(xx=FXX[index.GOSS,], yy=FY[index.GOSS,], 
                          beta, Var.a, Var.e, nc, R, p)
  GOSS.pred <- MSPE_fn(FY, FXX, FXX[index.GOSS,], FY[index.GOSS,], 
                       GOSS.Est[[5]], GOSS.Est[[6]], GOSS.Est[[7]], nc,C, R)
  GOSS.bt.mat <- GOSS.Est[[1]]
  GOSS.Var.a<- GOSS.Est[[2]]
  GOSS.Var.e<- GOSS.Est[[3]]
  GOSS.bt0.dif <- GOSS.Est[[4]]
  GOSS.bt <- GOSS.Est[[5]]
  ############################################################# estimate GIBOSS
  GIBOSS.Est <- Est_hat_cpp(xx=FXX[index.GIBOSS,], yy=FY[index.GIBOSS,], 
                            beta, Var.a, Var.e, nc.GIBOSS, R, p)
  GIBOSS.pred<- MSPE_fn(FY, FXX, FXX[index.GIBOSS,], FY[index.GIBOSS,], 
                        GIBOSS.Est[[5]], GIBOSS.Est[[6]], GIBOSS.Est[[7]], nc.GIBOSS,C, R)
  GIBOSS.bt.mat<- GIBOSS.Est[[1]]
  GIBOSS.Var.a<- GIBOSS.Est[[2]]
  GIBOSS.Var.e<- GIBOSS.Est[[3]]
  GIBOSS.bt0.dif<- GIBOSS.Est[[4]]
  GIBOSS.bt <- GIBOSS.Est[[5]]
  ############################################################# estimate CGOSS
  
  CGOSS.Est <- Est_hat_cpp(xx=FXX[final_index_CGOSS,], yy=FY[final_index_CGOSS,], 
                           beta, Var.a, Var.e, ncCGOSS, R_CGOSS, p)
  CGOSS.pred  <- MSPE_fn(FY.est, FX.est, FXX[final_index_CGOSS,], FY[final_index_CGOSS,], 
                         CGOSS.Est[[5]], CGOSS.Est[[6]], CGOSS.Est[[7]], ncCGOSS,C.est, R_CGOSS)
  CGOSS.bt.mat <- CGOSS.Est[[1]]
  CGOSS.Var.a<- CGOSS.Est[[2]]
  CGOSS.Var.e<- CGOSS.Est[[3]]
  CGOSS.bt0.dif <- CGOSS.Est[[4]]
  CGOSS.bt <- CGOSS.Est[[5]]
  ##########################################################  OSS
  EST_OSS_LM<-MSE_LM(FXX[index.OSS,],FY[index.OSS,],beta)
  OSS.pred<-OSS_mspe<-MSPE_LM(FXX,FY,EST_OSS_LM[[3]])
  OSS.bt.mat <- EST_OSS_LM[[2]]
  OSS.bt0.dif <- EST_OSS_LM[[1]]
  OSS.bt <- EST_OSS_LM[[3]]
  ############################################################# estimate IBOSS
  EST_IBOSS_LM<-MSE_LM(FXX[index.IBOSS,],FY[index.IBOSS,],beta)
  IBOSS.pred<-IBOSS_mspe<-MSPE_LM(FXX,FY,EST_IBOSS_LM[[3]])
  IBOSS.bt.mat <- EST_IBOSS_LM[[2]]
  IBOSS.bt0.dif <- EST_IBOSS_LM[[1]]
  IBOSS.bt <- EST_IBOSS_LM[[3]]
  
  #########IBOSSLMM
  IBOSS.Est.LMM <- Est_hat_cpp(xx=FXX[index.IBOSS,], yy=FY[index.IBOSS,], 
                               beta, Var.a, Var.e, nc3, R, p)
  IBOSS.pred.LMM <- MSPE_fn(FY, FXX, FXX[index.IBOSS,], FY[index.IBOSS,], 
                            IBOSS.Est.LMM[[5]], IBOSS.Est.LMM[[6]], IBOSS.Est.LMM[[7]], nc3, C, R)
  IBOSSLMM.bt.mat <- IBOSS.Est.LMM[[1]]
  IBOSSLMM.Var.a<- IBOSS.Est.LMM[[2]]
  IBOSSLMM.Var.e<- IBOSS.Est.LMM[[3]]
  IBOSSLMM.bt0.dif <- IBOSS.Est.LMM[[4]]
  IBOSSLMM.bt <- IBOSS.Est.LMM[[5]]
  
  #########OSSLMM
  OSS.Est.LMM <- Est_hat_cpp(xx=FXX[index.OSS,], yy=FY[index.OSS,], 
                             beta, Var.a, Var.e, nc2, R, p)
  OSS.pred.LMM <- MSPE_fn(FY, FXX, FXX[index.OSS,], FY[index.OSS,], 
                          OSS.Est.LMM[[5]], OSS.Est.LMM[[6]], OSS.Est.LMM[[7]], nc2, C, R)
  OSSLMM.bt.mat <- OSS.Est.LMM[[1]]
  OSSLMM.Var.a<- IBOSS.Est.LMM[[2]]
  OSSLMM.Var.e<- IBOSS.Est.LMM[[3]]
  OSSLMM.bt0.dif <- IBOSS.Est.LMM[[4]]
  OSSLMM.bt <- IBOSS.Est.LMM[[5]]
  
  
  ########################
  
  
  
  rec1<-cbind(CGOSS.bt.mat, IBOSS.bt.mat, OSS.bt.mat, GOSS.bt.mat,GIBOSS.bt.mat,IBOSSLMM.bt.mat,OSSLMM.bt.mat)
  rec2<-cbind(CGOSS.Var.a, GOSS.Var.a,GIBOSS.Var.a,IBOSSLMM.Var.a,OSSLMM.Var.a)
  rec3<-cbind(CGOSS.Var.e,GOSS.Var.e,GIBOSS.Var.e,IBOSSLMM.Var.e,OSSLMM.Var.e)
  rec4 <- cbind(CGOSS.pred, IBOSS.pred, OSS.pred, GOSS.pred,GIBOSS.pred,IBOSS.pred.LMM,OSS.pred.LMM)
  
  return(list(rec1,rec2,rec3,rec4))
}
#########################

## ==== CONFIG ====
file_path <- "climate_change_impact_on_agriculture_2024.csv"
y_col     <- "Crop_Yield_MT_per_HA"
group_col <- "Region"

## ==== LOAD ====
df <- read.csv(file_path, stringsAsFactors = FALSE)

## ==== CHECK ====
stopifnot(y_col %in% names(df), group_col %in% names(df))
df[[group_col]] <- factor(df[[group_col]])

## ==== Y ====
Y <- df[[y_col]]

## ==== X: numeric (NO scaling) + categorical (single-column numeric coding) ====
# numeric predictors (exclude Y)
num_cols <- names(df)[sapply(df, is.numeric)]
num_cols <- setdiff(num_cols, y_col)
X_num <- df[, num_cols, drop = FALSE]

# categorical predictors (exclude Y)
cat_cols <- names(df)[sapply(df, function(z) is.character(z) || is.factor(z))]
cat_cols <- setdiff(cat_cols, y_col)

# encode each categorical variable into ONE numeric column
encode_single_col <- function(x) {
  f <- factor(x)
  k <- nlevels(f)
  if (k == 1) {
    return(rep(0L, length(f)))
  } else if (k == 2) {
    return(ifelse(f == levels(f)[1], -1L, 1L))
  } else {
    # multi-class: symmetric integers centered at 0
    codes <- seq_len(k) - (k + 1) / 2
    return(as.integer(codes[as.integer(f)]))
  }
}

X_cat <- as.data.frame(
  lapply(df[, cat_cols, drop = FALSE], encode_single_col)
)

# final X (no scaling, no expansion)
X <- as.matrix(cbind(X_num, X_cat))
X <- apply(X, 2, as.numeric)

## ==== GROUPED PLACEMENT ====
# split by group
idx_by_group <- split(seq_len(nrow(df)), df[[group_col]])
X_by_group   <- lapply(idx_by_group, function(i) X[i, , drop = FALSE])
Y_by_group   <- lapply(idx_by_group, function(i) Y[i])

# sorted contiguous blocks + SCC
ord           <- order(df[[group_col]])
X_sorted      <- X[ord, , drop = FALSE]
Y_sorted      <- Y[ord]
group_sorted  <- df[[group_col]][ord]
group_sizes   <- as.integer(table(group_sorted))
SCC           <- c(0L, cumsum(group_sizes))



result = Comp(X=X_sorted,Y=Y_sorted,SSC=SCC,groupsize=group_sizes,n=1e3,setted_cluster=10)
result



