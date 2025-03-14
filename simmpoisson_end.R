##################
## Template for simulating multivariate GLMMs with Gaussian responses
##################
##################
## Tempkate for simulating multivaruate GLMMs with one random efffects
##################
rm(list = ls())
getwd()
#library(lme4)
library(scales)
library(ggplot2)
library(Matrix)
library(MCMCglmm)
library(mvtnorm) 
setwd("C:\\Users\\86188\\Desktop\\simulation for multi Poission response")

#setwd("D:\\OneDrive\\桌面\\composite_likelihood_and_variational_inference\\VI-PL\\VI-MGLMM")
#setwd("D:\\VI-MGLMM\\simulation for multi Poission response")

#source("aples-auxilaryfunctions.R")
source("CalculateB.Rs")
source("GVArandint.Rs")
source("mGLMM.Rs")

## Dataset generation for independent cluster GLMM; b_i ~ N(0,Sigma)
## Intercept must be manually included in X and Z if desigray
#这个函数的作用是生成多组观测数据，用于分析多级广义线性混合模型（MGLMM
gendat.mglmm <- function(id, fixed.MM, ran.MM, beta, Sigma, phi = NULL, cutoffs = NULL, family, save.MM = FALSE) {
  K <- nrow(beta)
  for(k in 1:K){
    #判断固定效应设计矩阵的第一列是否全部等于1，即判断是否包含截距项。
    if(!all(fixed.MM[[k]][,1] == 1)) 
      warning("Is your first column of fixed.MM an intercept? Should it be?")
    #检查随机效应设计矩阵ran.MM中每个响应变量对应的第一列是否为截距项（值全为1）
    if(!all(ran.MM[[k]][,1] == 1)) 
      warning("Is your first column of ran.MM an intercept? Should it be?")
  }
  #检查输入的family参数是否属于可允许的类型
  if(!(family %in% c("GAUSSIAN","POISSON","BINARY","ORDINAL"))) 
    stop("Specified family not permitted. Sorry!")
  for(k in 1:K){
    #语句检查fixed.MM[[k]]矩阵的列名是否为空。如果为空，
    #则通过paste函数生成列名，例如"X1"，"X2"等，并赋值给colnames(fixed.MM[[k]])
    if(is.null(colnames(fixed.MM[[k]]))) 
      colnames(fixed.MM[[k]]) <- paste("X", 1:ncol(fixed.MM[[k]]), sep = "")
    if(is.null(colnames(ran.MM[[k]]))) 
      colnames(ran.MM[[k]]) <- paste("Z", 1:ncol(ran.MM[[k]]), sep = "")
    #检查Sigma矩阵的行数是否与K*ncol(ran.MM[[k]])相等。如果不相等，则触发错误并显示错误信息
    #这是为了确保矩阵Sigma和ran.MM的维度一致。
    if(nrow(Sigma) != K*ncol(ran.MM[[k]])) 
      stop("Dimensions of Sigma, ran.MM are inconsistent. Please rectify. Thanks!")
    #检查beta矩阵的列数是否与fixed.MM[[k]]的列数相等
    if(ncol(beta) != ncol(fixed.MM[[k]])) 
      stop("Dimensions of beta and ran.MM are inconsistent. Please rectify. Thanks!")
  }
  m <- length(unique(id))
  randim <- 0
  for (i in 1:K){
    randim <- randim + ncol(ran.MM[[i]])
  }
  ##生成一个具有多元正态分布（以 Sigma+diag(x = 1e-8, nrow = nrow(Sigma)) 为协方差矩阵）的随机向量 true.b。
  ##其中，m 表示生成的随机向量的数量，rep(0,randim) 表示长度为 randim 的全零向量的均值向量
  true.b <- rmvnorm(m, rep(0,randim), sigma = Sigma+diag(x = 1e-8, nrow = nrow(Sigma)))
  true.b[,which(diag(Sigma) == 0)] <- 0
  ##创建一个名为 sim.y 的空矩阵，其行数为 fixed.MM[[1]] 的行数，列数为 K
  sim.y <- matrix(NA, nrow = nrow(fixed.MM[[1]]), ncol = K)
  for(i in 1:m) {
    sel.i <- which(id == i)
    for(k in 1:K){
      #计算当前个体的对应响应变量的预测值 eta，fixed.MM[[k]][sel.i,]%*%t(beta)[,k] 表示固定效应矩阵的乘积
      #ran.MM[[k]][sel.i,] 表示随机效应矩阵的对应部分
      #matrix(true.b[i,], nrow = K, ncol = ncol(ran.MM[[k]]), byrow = TRUE) 表示将随机效应向量转化为矩阵，用于与随机效应矩阵做乘积。
      #eta <- fixed.MM[[k]][sel.i,]%*%t(beta)[,k] + ran.MM[[k]][sel.i,]%*%t(matrix(true.b[i,], nrow = K, ncol = ncol(ran.MM[[k]]), byrow = TRUE))[,k] ## length(sel.i) by K 
      eta <- fixed.MM[[k]][sel.i,]%*%t(beta)[,k] + ran.MM[[k]][sel.i,]*t(matrix(true.b[i,], nrow = K, ncol = ncol(ran.MM[[k]]), byrow = TRUE))[,k] ## length(sel.i) by K 
      #对于高斯分布的响应变量，使用 rnorm() 函数生成随机数
      if(family == "GAUSSIAN") sim.y[sel.i,k] <- rnorm(length(sel.i), mean = eta, sd = sqrt(phi[k]))
      #对于泊松分布的响应变量，使用 rpois() 函数生成随机数
      pmean <- exp(eta)
      if(family == "POISSON" ) sim.y[sel.i,k] <- rpois(length(sel.i), pmean)
      #对于二项分布的响应变量，使用 rbinom() 函数生成随机数
      bmean <- 1/(1+exp(-eta))
      if(family == "BINARY"  ) sim.y[sel.i,k] <- rbinom(length(sel.i), 1, bmean) 
      if(family == "ORDINAL") sim.y[sel.i,k] <- rorddata(cutoffs = cutoffs[k,], eta = eta)
    }
  }
  
  #创建一个 NA 值填充的三维数组 b.mat，其尺寸为 (m, K, ncol(ran.MM[[1]]))
  b.mat <- array(NA, dim = c(m, K, ncol(ran.MM[[1]]))); 
  #使用 for 循环遍历所有个体，将该个体对应的随机效应向量转化为矩阵，并将其填充到 b.mat 中对应的位置上
  for(i in 1:m) b.mat[i,,] <- matrix(true.b[i,], nrow = K, ncol = ncol(ran.MM[[1]]), byrow = TRUE)	
  #创建一个 list 对象 out
  out <- list(y = sim.y, id = id, beta = beta, b = true.b, b.mat = b.mat, Sigma = Sigma, phi = phi, cutoffs = cutoffs, nonzero.fixef = (beta != 0), nonzero.ranef = matrix(diag(Sigma) != 0, nrow = K, byrow = TRUE))
  #如果需要保存固定效应和随机效应矩阵，将它们分别加入 out 中，并返回 out。如果不需要保存，直接返回 out
  if(save.MM) { out$fixed.MM <- fixed.MM; out$ran.MM <- ran.MM }
  return(out)
}


K <- 2  ## Noumber of responses
p <- 2 ## Noumber of covariates excluding intercept 
pr<- 1##随机效应矩阵的列数

true.beta <- matrix(c(-1,-0.5,0.5,0.3,-0.8,-1),ncol=3)
#true.beta <- matrix(c(-0.5,0.6,,-1,-0.5,0.8,-0.4),ncol=3)

true.Sigma <- matrix(c(1,0.5,0.5,0.8),ncol=2)# 真实的协方差矩阵 true.Sigma，对角线上的元素为 1，其他元素为 0.5
#true.Sigma <- matrix(c(0.8,0.7,0.7,0.5),ncol=2)
true.phi <- rep(1,K) ## True sigma2，真实的方差向量 true.phi，每个响应变量的方差都为 1,用于高斯分布生成数据
M <- c(10,50,100,500,1000,3000,5000,9000) ##subject
n <- 10 ##item
MS <- 100#100次模拟


beta01.gva  <- matrix(0,length(M)*MS,1); dim(beta01.gva)  <- c(length(M),1,MS)
beta11.gva  <- matrix(0,length(M)*MS,1); dim(beta11.gva)  <- c(length(M),1,MS)
beta12.gva  <- matrix(0,length(M)*MS,1); dim(beta12.gva)  <- c(length(M),1,MS)
beta02.gva  <- matrix(0,length(M)*MS,1); dim(beta02.gva)  <- c(length(M),1,MS)
beta21.gva  <- matrix(0,length(M)*MS,1); dim(beta21.gva)  <- c(length(M),1,MS)
beta22.gva  <- matrix(0,length(M)*MS,1); dim(beta22.gva)  <- c(length(M),1,MS)

sigma1.gva <- matrix(0,length(M)*MS,1); dim(sigma1.gva)  <- c(length(M),1,MS)
sigma2.gva <- matrix(0,length(M)*MS,1); dim(sigma2.gva)  <- c(length(M),1,MS)
rho.gva <- matrix(0,length(M)*MS,1); dim(rho.gva)  <- c(length(M),1,MS)

beta01.mcmc <- matrix(0,length(M)*MS,1); dim(beta01.mcmc)  <- c(length(M),1,MS)
beta11.mcmc <- matrix(0,length(M)*MS,1); dim(beta11.mcmc)  <- c(length(M),1,MS)
beta12.mcmc <- matrix(0,length(M)*MS,1); dim(beta12.mcmc)  <- c(length(M),1,MS)
beta02.mcmc <- matrix(0,length(M)*MS,1); dim(beta02.mcmc)  <- c(length(M),1,MS)
beta21.mcmc <- matrix(0,length(M)*MS,1); dim(beta21.mcmc)  <- c(length(M),1,MS)
beta22.mcmc <- matrix(0,length(M)*MS,1); dim(beta22.mcmc)  <- c(length(M),1,MS)

sigma1.mcmc <- matrix(0,length(M)*MS,1); dim(sigma1.mcmc)  <- c(length(M),1,MS)
sigma2.mcmc <- matrix(0,length(M)*MS,1); dim(sigma2.mcmc)  <- c(length(M),1,MS)
rho.mcmc <- matrix(0,length(M)*MS,1); dim(rho.mcmc)  <- c(length(M),1,MS)


times.gva <- matrix(0,length(M)*MS,3,1); dim(times.gva)  <- c(length(M),1,MS,3)
times.mcmc <- matrix(0,length(M)*MS,1); dim(times.mcmc)  <- c(length(M),1,MS)
converged.gva <- matrix(0,length(M)*MS,K); dim(converged.gva)  <- c(length(M),K,MS)
converged.mcmc <- matrix(0,length(M)*MS,1); dim(converged.mcmc)  <- c(length(M),1,MS)

ITER0 <- 1
JTER0 <- 1
trial0 <-  1

for (ITER in ITER0:length(M)) {
  ITER0  <- 1				
  for (JTER in JTER0:1) {
    JTER0  <- 1				
    for (trial in trial0:MS) {
      trial0 <- 1 
      set.seed(trial)
      p <- 2
      m <- M[ITER]
      id <- rep(1:m, each = n)
      true.Sigma <- matrix(c(1,0.2,0.2,1),ncol=2)
      
      
      fixed.MM <- list()
      ran.MM <- list()
      XSigma <- matrix(c(1,0,0,1),2,2)
      #生成一个大小为m*n的随机数矩阵X，这些随机数来自于均值为0、协方差矩阵为XSigma的多元正态分布。
      X1<- rmvnorm(m*n,rep(0,2),sigma=XSigma)
      X2<- rmvnorm(m*n,rep(0,2),sigma=XSigma)
      
      for(k in 1:K){
        fixed.MM[[k]] <- cbind(1,X1[,k],X2[,k])
        ran.MM[[k]] <- as.matrix(fixed.MM[[k]][,1],ncol=1)
      }
      beta <- true.beta
      Sigma <- true.Sigma
      phi = true.phi
      family = "POISSON"
      simy1 <- gendat.mglmm(id = id, fixed.MM = fixed.MM, ran.MM = ran.MM, beta = true.beta, Sigma = true.Sigma, phi = true.phi, family = "POISSON")
      
      ###Dec 31st,2022_start###
      ############################################
      vy <- matrix(simy1$y[,1],ncol=1)
      mX <- fixed.MM[[1]]																
      vn <- rep(n,m)		
      vbeta <- matrix(0,p+1,1)  
      sigma.true <- true.Sigma[1,1]#从true.Sigma矩阵中提取第一行第一列的元素，并将其赋值给sigma.true变量
      gamma <- log(2*sigma.true^2)  
      vmu   <- matrix(0,m,1)  
      vzeta <- matrix(1,m,1)
      
      ###############################################################################
      cat("trial=",trial,"    m=",m,"    sigma2=",sigma.true^2,"\n")
      cat("Fitting model using GVA \n")						
      
      bval1  <- proc.time()
      res <- try(GVA.fit1 <- GVA.randint.FIT(vbeta,exp(gamma),vmu,exp(vzeta),family,vy,mX,vn,id), silent = TRUE) 
      eval1  <- proc.time()
      if((class(res)[1])!="try-error"){ 
        
        beta01.gva[ITER,JTER,trial] <- res$vbeta[1]
        beta11.gva[ITER,JTER,trial] <- res$vbeta[2]
        beta12.gva[ITER,JTER,trial] <- res$vbeta[3]
        sigma1.gva[ITER,JTER,trial] <- res$sigma2
        
        times.gva[ITER,JTER,trial,1] <- eval1[3] - bval1[3]
        converged.gva[ITER,JTER,trial] <- TRUE
      } else {
        converged.gva[ITER,JTER,trial] <- FALSE#表明拟合未成功
      }
      
      vy <- matrix(simy1$y[,2],ncol=1)
      mX <- fixed.MM[[2]]															
      vbeta <- matrix(0,p+1,1)  
      sigma.true <- true.Sigma[2,2]
      gamma <- log(2*sigma.true^2)  
      vmu   <- matrix(0,m,1)  
      vzeta <- matrix(1,m,1)
      ###############################################################################
      cat("trial=",trial,"    m=",m,"    sigma2=",sigma.true,"\n")
      cat("Fitting model using GVA \n")						
      bval2  <- proc.time()
      res <- try(GVA.fit2 <- GVA.randint.FIT(vbeta,exp(gamma),vmu,exp(vzeta),family,vy,mX,vn,id), silent = TRUE) 
      eval2  <- proc.time()
      
      if((class(res)[1])!="try-error"){ 
        
        beta02.gva[ITER,JTER,trial] <- res$vbeta[1]
        beta21.gva[ITER,JTER,trial] <- res$vbeta[2]
        beta22.gva[ITER,JTER,trial] <- res$vbeta[3]
        sigma2.gva[ITER,JTER,trial] <- res$sigma2
        
        times.gva[ITER,JTER,trial,2] <- eval2[3] - bval2[3]
        converged.gva[ITER,JTER,trial] <- TRUE
        
      } else {
        converged.gva[ITER,JTER,trial] <- FALSE
      }
      
      bval3  <- proc.time()
      vcov <- cov(simy1$y)[1,2]
      vde <- exp((sigma1.gva[ITER,JTER,trial]+sigma2.gva[ITER,JTER,trial])/2)*mean(exp((fixed.MM[[1]]%*%c(beta01.gva[ITER,JTER,trial],beta11.gva[ITER,JTER,trial],beta12.gva[ITER,JTER,trial])+fixed.MM[[2]]%*%c(beta02.gva[ITER,JTER,trial],beta21.gva[ITER,JTER,trial],beta22.gva[ITER,JTER,trial]))))
      rho.gva[ITER,JTER,trial] <- log(1+vcov/vde)/(sqrt(sigma1.gva[ITER,JTER,trial])*sqrt(sigma2.gva[ITER,JTER,trial]))
      eval3  <- proc.time()
      times.gva[ITER,JTER,trial,3] <- eval3[3] - bval3[3]
      
      
      dat1 <- data.frame(
        vy1 = c(simy1$y[,1]),
        mX11 = c(fixed.MM[[1]][,1]),
        mX12 = c(fixed.MM[[1]][,2]),
        mX13 = c(fixed.MM[[1]][,3]),
        factor = factor(id),
        ID = c(1:length(id))
      )
      
      dat2 <- data.frame(
        vy2 = c(simy1$y[,2]),
        mX11 = c(fixed.MM[[2]][,1]),
        mX12 = c(fixed.MM[[2]][,2]),
        mX13 = c(fixed.MM[[2]][,3]),
        factor = factor(id),
        ID = c(1:length(id))
      )
      
      dat2$ID <- dim(dat1)[1] + 1:length(id)
      DAT <- merge(dat1,dat2,all = TRUE)
      #prior = list(R = list(V = diag(2), nu = 0, fix = 2), G = list(G1 = list(V = diag(2), nu = 1)))
      prior = list(
        R = list(V = diag(2), nu = 0.002),  # 较弱的先验
        G = list(G1 = list(V = diag(2), nu = 0.002))
      )
      bval4  <- proc.time()
      #res <- try(MCMCmGLMM.fit <- MCMCglmm(cbind(vy1,vy2) ~trait-1 + trait:mX12+trait:mX13, random = ~us(trait):factor,
                                          # rcov = ~idh(trait):units, family = c("poisson", "poisson"),
                                           #data = DAT, prior = prior, verbose = FALSE), silent = TRUE) 
      res <- try(MCMCglmm(cbind(vy1, vy2) ~ trait - 1 + trait:mX12 + trait:mX13, 
                          random = ~us(trait):factor, 
                          rcov = ~idh(trait):units, 
                          family = c("poisson", "poisson"), 
                          data = DAT, 
                          prior = prior, 
                          nitt = 13000, 
                          burnin = 3000, 
                          thin = 10, 
                          verbose = FALSE), 
                 silent = TRUE)
      eval4  <- proc.time()
      
      if((class(res)[1])!="try-error"){ 
        beta01.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$solutions[1,1])
        beta11.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$solutions[3,1])
        beta12.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$solutions[5,1])
        beta02.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$solutions[2,1])
        beta21.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$solutions[4,1])
        beta22.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$solutions[6,1])
        
        sigma1.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$Gcovariances[1,1])
        sigma2.mcmc[ITER,JTER,trial] <- as.numeric(summary(res)$Gcovariances[4,1])
        rho.mcmc[ITER,JTER,trial]   <- as.numeric(summary(res)$Gcovariances[2,1])/sqrt(as.numeric(summary(res)$Gcovariances[1,1])*as.numeric(summary(res)$Gcovariances[4,1]))
        
        times.mcmc[ITER,JTER,trial] <- eval4[3] - bval4[3]
      }
    }
  }
}

meantime.gva <- c()
meantime.mcmc <- c()
meanbeta01.gva <- c()
meanbeta01.mcmc <- c()
meanbeta11.gva <- c()
meanbeta11.mcmc <- c()
meanbeta12.gva <- c()
meanbeta12.mcmc <-c()
meanbeta02.gva <- c()
meanbeta02.mcmc <- c()
meanbeta21.gva <- c()
meanbeta21.mcmc <- c()
meanbeta22.gva <- c()
meanbeta22.mcmc <- c()
meansigma1.gva <- c()
meansigma1.mcmc <- c()
meansigma2.gva <- c()
meansigma2.mcmc <- c()
meanrho.gva <- c()
meanrho.mcmc <- c()
mse.beta01.gva <-c()
mse.beta11.gva <-c()
mse.beta12.gva <-c()

mse.beta02.gva <-c()
mse.beta21.gva <-c()
mse.beta22.gva <-c()

mse.sigma1.gva <-c()
mse.sigma2.gva <-c()
mse.rho.gva <-c()

mse.beta01.mcmc <-c()
mse.beta11.mcmc <-c()
mse.beta12.mcmc <-c()

mse.beta02.mcmc <-c()
mse.beta21.mcmc <-c()
mse.beta22.mcmc <-c()

mse.sigma1.mcmc <-c()
mse.sigma2.mcmc <-c()
mse.rho.mcmc <-c()
indOutliers <- function(x,nsd=3) {
  z.scores <- (x - median(x))/mad(x)
  return(which(abs(z.scores)>nsd))
}
for (ITER in 1:length(M)) {
  for (JTER in 1:1) {
    
    m <- M[ITER]
    index0 <- which(converged.gva[ITER,JTER,] & converged.gva[ITER,JTER,]) 
    all_outliers <- list(
      indOutliers(beta01.gva[ITER,JTER,index0]),
      indOutliers(beta11.gva[ITER,JTER,index0]),
      indOutliers(beta12.gva[ITER,JTER,index0]),
      indOutliers(beta02.gva[ITER,JTER,index0]),
      indOutliers(beta21.gva[ITER,JTER,index0]),
      indOutliers(beta22.gva[ITER,JTER,index0]),
      indOutliers(beta01.mcmc[ITER,JTER,index0]),
      indOutliers(beta11.mcmc[ITER,JTER,index0]),
      indOutliers(beta12.mcmc[ITER,JTER,index0]),
      indOutliers(beta02.mcmc[ITER,JTER,index0]),
      indOutliers(beta21.mcmc[ITER,JTER,index0]),
      indOutliers(beta22.mcmc[ITER,JTER,index0]),
      indOutliers(sigma1.gva[ITER,JTER,index0]),
      indOutliers(sigma2.gva[ITER,JTER,index0]),
      indOutliers(sigma1.mcmc[ITER,JTER,index0]),
      indOutliers(sigma2.mcmc[ITER,JTER,index0]),
      indOutliers(rho.gva[ITER,JTER,index0]),
      indOutliers(rho.mcmc[ITER,JTER,index0])
    )
    index0.outliers <- unique(unlist(all_outliers))               
    if (length(index0.outliers)>0) {              
      index0 <- index0[-index0.outliers]     
    }
    
    print(length(index0))
  
    dp <- 5
    # 计算均值
    meantime_gva <- mean(times.gva[ITER,JTER,index0,1] + times.gva[ITER,JTER,index0,2] + times.gva[ITER,JTER,index0,3])
    meantime_mcmc <- mean(times.mcmc[ITER,JTER,index0])
    meantime.gva[[ITER]] <- meantime_gva
    meantime.mcmc[[ITER]] <- meantime_mcmc
    # 存储GVA均值
    meanbeta01.gva[[ITER]] <- round(mean(beta01.gva[ITER,JTER,index0]), dp)
    meanbeta11.gva[[ITER]] <- round(mean(beta11.gva[ITER,JTER,index0]), dp)
    meanbeta12.gva[[ITER]] <- round(mean(beta12.gva[ITER,JTER,index0]), dp)
    meanbeta02.gva[[ITER]] <- round(mean(beta02.gva[ITER,JTER,index0]), dp)
    meanbeta21.gva[[ITER]] <- round(mean(beta21.gva[ITER,JTER,index0]), dp)
    meanbeta22.gva[[ITER]] <- round(mean(beta22.gva[ITER,JTER,index0]), dp)
    meansigma1.gva[[ITER]] <- round(mean(sigma1.gva[ITER,JTER,index0]), dp)
    meansigma2.gva[[ITER]] <- round(mean(sigma2.gva[ITER,JTER,index0]), dp)
    meanrho.gva[[ITER]] <- round(mean(rho.gva[ITER,JTER,index0]), dp)
    #存储GVAmse
    mse.beta01.gva[[ITER]]  <- round(mean((beta01.gva[ITER,JTER,index0] - beta[1,1])^2),dp)
    mse.beta11.gva[[ITER]]  <- round(mean((beta11.gva[ITER,JTER,index0] - beta[1,2])^2),dp)
    mse.beta12.gva[[ITER]]  <- round(mean((beta12.gva[ITER,JTER,index0] - beta[1,3])^2),dp)
    
    mse.beta02.gva[[ITER]] <- round(mean((beta02.gva[ITER,JTER,index0] - beta[2,1])^2),dp)
    mse.beta21.gva[[ITER]] <- round(mean((beta21.gva[ITER,JTER,index0] - beta[2,2])^2),dp)
    mse.beta22.gva[[ITER]] <- round(mean((beta22.gva[ITER,JTER,index0] - beta[2,3])^2),dp)
    
    mse.sigma1.gva[[ITER]]  <- round(mean((sigma1.gva[ITER,JTER,index0] - true.Sigma[1,1])^2),dp)
    mse.sigma2.gva[[ITER]]  <- round(mean((sigma2.gva[ITER,JTER,index0] - true.Sigma[2,2])^2),dp)
    mse.rho.gva[[ITER]] <- round(mean((rho.gva[ITER,JTER,index0] - (true.Sigma[1,2]/(sqrt(true.Sigma[1,1])*sqrt(true.Sigma[2,2]))))^2),dp)
    
    # 存储MCMC均值
    meanbeta01.mcmc[[ITER]] <- round(mean(beta01.mcmc[ITER,JTER,index0]), dp)
    meanbeta11.mcmc[[ITER]] <- round(mean(beta11.mcmc[ITER,JTER,index0]), dp)
    meanbeta12.mcmc[[ITER]] <- round(mean(beta12.mcmc[ITER,JTER,index0]), dp)
    meanbeta02.mcmc[[ITER]] <- round(mean(beta02.mcmc[ITER,JTER,index0]), dp)
    meanbeta21.mcmc[[ITER]] <- round(mean(beta21.mcmc[ITER,JTER,index0]), dp)
    meanbeta22.mcmc[[ITER]] <- round(mean(beta22.mcmc[ITER,JTER,index0]), dp)
    meansigma1.mcmc[[ITER]] <- round(mean(sigma1.mcmc[ITER,JTER,index0]), dp)
    meansigma2.mcmc[[ITER]] <- round(mean(sigma2.mcmc[ITER,JTER,index0]), dp)
    meanrho.mcmc[[ITER]] <- round(mean(rho.mcmc[ITER,JTER,index0]), dp)
    #存储MCMCmse
    mse.beta01.mcmc[[ITER]]  <- round(mean((beta01.mcmc[ITER,JTER,index0] - beta[1,1])^2),dp)
    mse.beta11.mcmc[[ITER]]  <- round(mean((beta11.mcmc[ITER,JTER,index0] - beta[1,2])^2),dp)
    mse.beta12.mcmc[[ITER]]  <- round(mean((beta12.mcmc[ITER,JTER,index0] - beta[1,3])^2),dp)
    
    mse.beta02.mcmc[[ITER]] <- round(mean((beta02.mcmc[ITER,JTER,index0] - beta[2,1])^2),dp)
    mse.beta21.mcmc[[ITER]] <- round(mean((beta21.mcmc[ITER,JTER,index0] - beta[2,2])^2),dp)
    mse.beta22.mcmc[[ITER]] <- round(mean((beta22.mcmc[ITER,JTER,index0] - beta[2,3])^2),dp)
    
    mse.sigma1.mcmc[[ITER]]  <- round(mean((sigma1.mcmc[ITER,JTER,index0] - true.Sigma[1,1])^2),dp)
    mse.sigma2.mcmc[[ITER]]  <- round(mean((sigma2.mcmc[ITER,JTER,index0] - true.Sigma[2,2])^2),dp)
    mse.rho.mcmc[[ITER]] <- round(mean((rho.mcmc[ITER,JTER,index0] - (true.Sigma[1,2]/(sqrt(true.Sigma[1,1])*sqrt(true.Sigma[2,2]))))^2),dp)
    
  }
}
###GVA-split###
print(meanbeta01.gva)
print(meanbeta11.gva)
print(meanbeta12.gva)
print(meanbeta02.gva)
print(meanbeta21.gva)
print(meanbeta22.gva)
print(meansigma1.gva)
print(meansigma2.gva)
print(meanrho.gva)
###MCMC###
print(meanbeta01.mcmc)
print(meanbeta11.mcmc)
print(meanbeta12.mcmc)
print(meanbeta02.mcmc)
print(meanbeta21.mcmc)
print(meanbeta22.mcmc)
print(meansigma1.mcmc)
print(meansigma2.mcmc)
print(meanrho.mcmc)
#GVA_mse
print(mse.beta01.gva)
print(mse.beta11.gva)
print(mse.beta12.gva)
print(mse.beta02.gva)
print(mse.beta21.gva)
print(mse.beta22.gva)
print(mse.sigma1.gva)
print(mse.sigma2.gva)
print(mse.rho.gva)
#MCMC_mse
print(mse.beta01.mcmc)
print(mse.beta11.mcmc)
print(mse.beta12.mcmc)
print(mse.beta02.mcmc)
print(mse.beta21.mcmc)
print(mse.beta22.mcmc)
print(mse.sigma1.mcmc)
print(mse.sigma2.mcmc)
print(mse.rho.mcmc)
data11<- data.frame(x=M*n,y=unlist(meantime.gva))
data12<- data.frame(x=M*n,y=unlist(meantime.mcmc))

data21<- data.frame(x=M*n,y=unlist(mse.beta01.gva))
data22<- data.frame(x=M*n,y=unlist(mse.beta01.mcmc))
data31<- data.frame(x=M*n,y=unlist(mse.beta11.gva))
data32<- data.frame(x=M*n,y=unlist(mse.beta11.mcmc))
data41<- data.frame(x=M*n,y=unlist(mse.beta12.gva))
data42<- data.frame(x=M*n,y=unlist(mse.beta12.mcmc))
data51<- data.frame(x=M*n,y=unlist(mse.beta02.gva))
data52<- data.frame(x=M*n,y=unlist(mse.beta02.mcmc))
data61<- data.frame(x=M*n,y=unlist(mse.beta21.gva))
data62<- data.frame(x=M*n,y=unlist(mse.beta21.mcmc))
data71<- data.frame(x=M*n,y=unlist(mse.beta22.gva))
data72<- data.frame(x=M*n,y=unlist(mse.beta22.mcmc))
data81<- data.frame(x=M*n,y=unlist(mse.sigma1.gva))
data82<- data.frame(x=M*n,y=unlist(mse.sigma1.mcmc))
data91<- data.frame(x=M*n,y=unlist(mse.sigma2.gva))
data92<- data.frame(x=M*n,y=unlist(mse.sigma2.mcmc))
data10_1<- data.frame(x=M*n,y=unlist(mse.rho.gva))
data10_2<- data.frame(x=M*n,y=unlist(mse.rho.mcmc))
# 自定义标签函数
format_labels <- function(x) {
  sapply(x, function(val) {
    if (val == 0) return("0")
    power <- floor(log10(val))
    base <- val / 10^power
    paste0(base, " * 10^", power)
  })
}

# 对x和y取对数
data11$log_x <- log10(data11$x)
data11$log_y <- log10(data11$y)
data12$log_x <- log10(data12$x)
data12$log_y <- log10(data12$y)
lm11<-lm(log_y ~ log_x, data = data11)
slope11<-summary(lm11)$coefficients[2,1]
lm12<-lm(log_y ~ log_x, data = data12)
slope12<-summary(lm12)$coefficients[2,1]
plot_time <- ggplot() +
  geom_point(data = data11, aes(x = log10(x), y = log10(y), color = "blue"), size = 3) +
  geom_point(data = data12, aes(x = log10(x), y = log10(y), color = "red"), size = 3) +
  geom_smooth(data = data11, aes(x = log10(x), y = log10(y), color = "blue"), method = "lm", se = FALSE) +
  geom_smooth(data = data12, aes(x = log10(x), y = log10(y), color = "red"), method = "lm", se = FALSE) +
  annotate("text", x = log10(max(data11$x)) * 0.8, y = log10(max(data11$y)) * 0.8, 
           label = paste("Slope1:", round(coef(lm11)[2], 5)), color = "blue", size = 4) +
  annotate("text", x = log10(max(data12$x)) * 0.8, y = log10(max(data12$y)) * 0.9, 
           label = paste("Slope2:", round(coef(lm12)[2], 5)), color = "red", size = 4) +
  scale_color_manual(values = c("blue", "red"), labels = c("GVA", "MCMC")) +
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title = "Computational Cost", x = "Sample Size", y = "Log(Time/s)") +
  theme_minimal() +
  theme(
    legend.position = "top"
  )

data21$log_x <- log10(data21$x)
data21$log_y <- log10(data21$y)
data22$log_x <- log10(data22$x)
data22$log_y <- log10(data22$y)
lm21<-lm(log_y~log_x,data=data21)
slope21<-summary(lm21)$coefficients[2,1]
lm22<-lm(log_y~log_x,data=data22)
slope22<-summary(lm22)$coefficients[2,1]
plot_beta01<-ggplot()+
  geom_point(data = data21,aes(x = log10(x),y = log10(y),color ="blue"), size =3)+	
  geom_point(data = data22,aes(x = log10(x),y = log10(y),color="red"), size =3)+	
  geom_smooth(data = data21,aes(x = log10(x),y = log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data =data22,aes(x = log10(x),y= log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x = log10(3000),y= log10(0.02), label = paste("Slope1:",round(slope21,4)),color ="blue",size=4)+
  annotate("text",x = log10(3000),y=log10(0.04), label = paste("Slope2:",round(slope22,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size ", y=expression(Log("MSE of "* hat(beta)[0*1])))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data31$log_x <- log10(data31$x)
data31$log_y <- log10(data31$y)
data32$log_x <- log10(data32$x)
data32$log_y <- log10(data32$y)
lm31<-lm(log_y~log_x,data=data31)
slope31<-summary(lm31)$coefficients[2,1]
lm32<-lm(log_y~log_x,data=data32)
slope32<-summary(lm32)$coefficients[2,1]
plot_beta11<-ggplot()+
  geom_point(data = data31,aes(x = log10(x), y = log10(y), color = "blue"), size = 3) +
  geom_point(data = data32,aes(x = log10(x), y = log10(y), color = "red"), size = 3) +
  geom_smooth(data = data31,aes(x = log10(x), y = log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data = data32,aes(x = log10(x), y = log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x = log10(3000),y= log10(0.0035), label = paste("Slope1:",round(slope31,4)),color ="blue",size=4)+
  annotate("text",x = log10(3000),y= log10(0.006), label = paste("Slope2:",round(slope32,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size", y=expression(Log("MSE of "* hat(beta[11]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data41$log_x <- log10(data41$x)
data41$log_y <- log10(data41$y)
data42$log_x <- log10(data42$x)
data42$log_y <- log10(data42$y)
lm41<-lm(log_y~log_x,data=data41)
slope41<-summary(lm41)$coefficients[2,1]
lm42<-lm(log_y~log_x,data=data42)
slope42<-summary(lm42)$coefficients[2,1]
plot_beta12<-ggplot()+
  geom_point(data = data41,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data42,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data41,aes(x = log10(x),y = log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data = data42,aes(x = log10(x),y = log10(y),color ="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(8000),y= log10(0.0005), label = paste("Slope1:",round(slope41,4)),color ="blue",size=4)+
  annotate("text",x =log10(8000),y=log10(0.0025), label = paste("Slope2:",round(slope42,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size", y=expression(Log("MSE of "* hat(beta[12]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data51$log_x <- log10(data51$x)
data51$log_y <- log10(data51$y)
data52$log_x <- log10(data52$x)
data52$log_y <- log10(data52$y)
lm51<-lm(log_y~log_x,data=data51)
slope51<-summary(lm51)$coefficients[2,1]
lm52<-lm(log_y~log_x,data=data52)
slope52<-summary(lm52)$coefficients[2,1]
plot_beta02<-ggplot()+
  geom_point(data = data51,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data52,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data51,aes(x=log10(x),y=log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data = data52,aes(x=log10(x),y=log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(2000),y= log10(0.045), label = paste("Slope1:",round(slope51,4)),color ="blue",size=4)+
  annotate("text",x =log10(2000),y=log10(0.08), label = paste("Slope2:",round(slope52,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size", y=expression(Log("MSE of "*hat(beta[0*2]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data61$log_x <- log10(data61$x)
data61$log_y <- log10(data61$y)
data62$log_x <- log10(data62$x)
data62$log_y <- log10(data62$y)
lm61<-lm(log_y~log_x,data=data61)
slope61<-summary(lm61)$coefficients[2,1]
lm62<-lm(log_y~log_x,data=data62)
slope62<-summary(lm62)$coefficients[2,1]
plot_beta21<-ggplot()+
  geom_point(data = data61,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data62,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data61,aes(x =log10(x),y =log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data =data62,aes(x=log10(x),y=log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(3000),y=log10(0.0025), label = paste("Slope1:",round(slope61,4)),color ="blue",size=4)+
  annotate("text",x =log10(3000),y=log10(0.005), label = paste("Slope2:",round(slope62,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size ", y=expression(Log("MSE of "*hat(beta[21]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data71$log_x <- log10(data71$x)
data71$log_y <- log10(data71$y)
data72$log_x <- log10(data72$x)
data72$log_y <- log10(data72$y)
lm71<-lm(log_y~log_x,data=data71)
slope71<-summary(lm71)$coefficients[2,1]
lm72<-lm(log_y~log_x,data=data72)
slope72<-summary(lm72)$coefficients[2,1]
plot_beta22<-ggplot()+
  geom_point(data = data71,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data72,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data71,aes(x =log10(x),y =log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data =data72,aes(x=log10(x),y= log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(3000),y= log10(0.001), label = paste("Slope1:",round(slope71,4)),color ="blue",size=4)+
  annotate("text",x =log10(3000),y=log10(0.004), label = paste("Slope2:",round(slope72,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size", y=expression(Log("MSE of "*hat(beta[22]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data81$log_x <- log10(data81$x)
data81$log_y <- log10(data81$y)
data82$log_x <- log10(data82$x)
data82$log_y <- log10(data82$y)
lm81<-lm(log_y~log_x,data=data81)
slope81<-summary(lm81)$coefficients[2,1]
lm82<-lm(log_y~log_x,data=data82)
slope82<-summary(lm82)$coefficients[2,1]
plot_sigma1<-ggplot()+
  geom_point(data = data81,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data82,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data81,aes(x =log10(x),y =log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data =data82,aes(x=log10(x),y= log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(1000),y= log10(0.002), label = paste("Slope1:",round(slope81,4)),color ="blue",size=4)+
  annotate("text",x =log10(10000),y=log10(0.1), label = paste("Slope2:",round(slope82,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size", y=expression(Log("MSE of "*hat(sigma[1]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data91$log_x <- log10(data91$x)
data91$log_y <- log10(data91$y)
data92$log_x <- log10(data92$x)
data92$log_y <- log10(data92$y)
lm91<-lm(log_y~log_x,data=data91)
slope91<-summary(lm91)$coefficients[2,1]
lm92<-lm(log_y~log_x,data=data92)
slope92<-summary(lm92)$coefficients[2,1]
plot_sigma2<-ggplot()+
  geom_point(data = data91,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data92,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data91,aes(x =log10(x),y =log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data =data92,aes(x=log10(x),y=log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(10000),y= log10(0.05), label = paste("Slope1:",round(slope91,4)),color ="blue",size=4)+
  annotate("text",x =log10(10000),y=log10(0.15), label = paste("Slope2:",round(slope92,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size", y=expression(Log("MSE of "*hat(sigma[2]))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

data10_1$log_x <- log10(data10_1$x)
data10_1$log_y <- log10(data10_1$y)
data10_2$log_x <- log10(data10_2$x)
data10_2$log_y <- log10(data10_2$y)
lm10_1<-lm(log_y~log_x,data=data10_1)
slope10_1<-summary(lm10_1)$coefficients[2,1]
lm10_2<-lm(log_y~log_x,data=data10_2)
slope10_2<-summary(lm10_2)$coefficients[2,1]
plot_rho<-ggplot()+
  geom_point(data = data10_1,aes(x =log10(x),y=log10(y),color ="blue"), size =3)+	
  geom_point(data = data10_2,aes(x=log10(x),y=log10(y),color="red"), size =3)+	
  geom_smooth(data = data10_1,aes(x =log10(x),y =log10(y),color ="blue"), method = "lm", se=FALSE)+	
  geom_smooth(data =data10_2,aes(x=log10(x),y=log10(y),color="red"),method ="lm",se=FALSE)+	
  annotate("text",x =log10(10000),y= log10(0.1), label = paste("Slope1:",round(slope10_1,4)),color ="blue",size=4)+
  annotate("text",x =log10(3000),y=log10(0.001), label = paste("Slope2:",round(slope10_2,4)),color = "red",size=4)+
  scale_color_manual(values = c("blue","red"),labels = c("GVA","MCMC"))+	
  scale_x_continuous(
    breaks = log10(c(1e2, 1e3, 1e4, 1e5)), 
    labels = scales::label_math(10^.x)
  ) + 
  labs(title="Mean Square Error",x="Sample size ", y=expression(Log("MSE of "*hat(rho))))+	
  theme_minimal()+	
  theme(
    legend.position = "top"
  )

for (ITER in 1:length(M)) {
  for (JTER in 1:1) {
    
    m <- M[ITER]
    index0 <- which(converged.gva[ITER,JTER,] & converged.gva[ITER,JTER,]) 
    all_outliers <- list(
      indOutliers(beta01.gva[ITER,JTER,index0]),
      indOutliers(beta11.gva[ITER,JTER,index0]),
      indOutliers(beta12.gva[ITER,JTER,index0]),
      indOutliers(beta02.gva[ITER,JTER,index0]),
      indOutliers(beta21.gva[ITER,JTER,index0]),
      indOutliers(beta22.gva[ITER,JTER,index0]),
      indOutliers(beta01.mcmc[ITER,JTER,index0]),
      indOutliers(beta11.mcmc[ITER,JTER,index0]),
      indOutliers(beta12.mcmc[ITER,JTER,index0]),
      indOutliers(beta02.mcmc[ITER,JTER,index0]),
      indOutliers(beta21.mcmc[ITER,JTER,index0]),
      indOutliers(beta22.mcmc[ITER,JTER,index0]),
      indOutliers(sigma1.gva[ITER,JTER,index0]),
      indOutliers(sigma2.gva[ITER,JTER,index0]),
      indOutliers(sigma1.mcmc[ITER,JTER,index0]),
      indOutliers(sigma2.mcmc[ITER,JTER,index0]),
      indOutliers(rho.gva[ITER,JTER,index0]),
      indOutliers(rho.mcmc[ITER,JTER,index0])
    )
    index0.outliers <- unique(unlist(all_outliers))               
    if (length(index0.outliers)>0) {              
      index0 <- index0[-index0.outliers]     
    }
  }
}
data_gva_beta01 <- data.frame(
  group = c(
    rep("1*10^2", length(beta01.gva[1, , index0])),
    rep("5*10^2", length(beta01.gva[2, , index0])),
    rep("1*10^3", length(beta01.gva[3, , index0])),
    rep("5*10^3", length(beta01.gva[4, , index0])),
    rep("1*10^4", length(beta01.gva[5, , index0])),
    rep("3*10^4", length(beta01.gva[6, , index0])),
    rep("5*10^4", length(beta01.gva[7, , index0])),
    rep("9*10^4", length(beta01.gva[8, , index0]))
  ),
  value = c(
    beta01.gva[1, , index0],
    beta01.gva[2, , index0],
    beta01.gva[3, , index0],
    beta01.gva[4, , index0],
    beta01.gva[5, , index0],
    beta01.gva[6, , index0],
    beta01.gva[7, , index0],
    beta01.gva[8, , index0]
  )
)
data_gva_beta01$group <- factor(data_gva_beta01$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta01_gva <- ggplot(data_gva_beta01, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray", outlier.color = "gray", outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8) +
  labs(x = "Sample Size ", y = expression(hat(beta[0*1])), title = expression("GVA: Boxplots of " * hat(beta)[0 * 1])) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-2, 0)  # 限定纵坐标范围

data_gva_beta11 <- data.frame(
  group = c(
    rep("1*10^2", length(beta11.gva[1, , index0])),
    rep("5*10^2", length(beta11.gva[2, , index0])),
    rep("1*10^3", length(beta11.gva[3, , index0])),
    rep("5*10^3", length(beta11.gva[4, , index0])),
    rep("1*10^4", length(beta11.gva[5, , index0])),
    rep("3*10^4", length(beta11.gva[6, , index0])),
    rep("5*10^4", length(beta11.gva[7, , index0])),
    rep("9*10^4", length(beta11.gva[8, , index0]))
  ),
  value = c(
    beta11.gva[1, , index0],
    beta11.gva[2, , index0],
    beta11.gva[3, , index0],
    beta11.gva[4, , index0],
    beta11.gva[5, , index0],
    beta11.gva[6, , index0],
    beta11.gva[7, , index0],
    beta11.gva[8, , index0]
  )
)
data_gva_beta11$group <- factor(data_gva_beta11$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta11_gva <- ggplot(data_gva_beta11, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[11])), title = expression("GVA: Boxplots of " * hat(beta)[11]))+
  theme_minimal() +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.5, -0.5) 
data_gva_beta12 <- data.frame(
  group = c(
    rep("1*10^2", length(beta12.gva[1, , index0])),
    rep("5*10^2", length(beta12.gva[2, , index0])),
    rep("1*10^3", length(beta12.gva[3, , index0])),
    rep("5*10^3", length(beta12.gva[4, , index0])),
    rep("1*10^4", length(beta12.gva[5, , index0])),
    rep("3*10^4", length(beta12.gva[6, , index0])),
    rep("5*10^4", length(beta12.gva[7, , index0])),
    rep("9*10^4", length(beta12.gva[8, , index0]))
  ),
  value = c(
    beta12.gva[1, , index0],
    beta12.gva[2, , index0],
    beta12.gva[3, , index0],
    beta12.gva[4, , index0],
    beta12.gva[5, , index0],
    beta12.gva[6, , index0],
    beta12.gva[7, , index0],
    beta12.gva[8, , index0]
  )
)
data_gva_beta12$group <- factor(data_gva_beta12$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta12_gva <- ggplot(data_gva_beta12, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[12])), title = expression("GVA: Boxplots of " * hat(beta)[12]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.25, -0.75) 
data_gva_beta02 <- data.frame(
  group = c(
    rep("1*10^2", length(beta02.gva[1, , index0])),
    rep("5*10^2", length(beta02.gva[2, , index0])),
    rep("1*10^3", length(beta02.gva[3, , index0])),
    rep("5*10^3", length(beta02.gva[4, , index0])),
    rep("1*10^4", length(beta02.gva[5, , index0])),
    rep("3*10^4", length(beta02.gva[6, , index0])),
    rep("5*10^4", length(beta02.gva[7, , index0])),
    rep("9*10^4", length(beta02.gva[8, , index0]))
  ),
  value = c(
    beta02.gva[1, , index0],
    beta02.gva[2, , index0],
    beta02.gva[3, , index0],
    beta02.gva[4, , index0],
    beta02.gva[5, , index0],
    beta02.gva[6, , index0],
    beta02.gva[7, , index0],
    beta02.gva[8, , index0]
  )
)
data_gva_beta02$group <- factor(data_gva_beta02$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta02_gva <- ggplot(data_gva_beta02, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[0*2])), title = expression("GVA: Boxplots of " * hat(beta)[0*2]))+
  theme_minimal() +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-2, 0) 
data_gva_beta21 <- data.frame(
  group = c(
    rep("1*10^2", length(beta21.gva[1, , index0])),
    rep("5*10^2", length(beta21.gva[2, , index0])),
    rep("1*10^3", length(beta21.gva[3, , index0])),
    rep("5*10^3", length(beta21.gva[4, , index0])),
    rep("1*10^4", length(beta21.gva[5, , index0])),
    rep("3*10^4", length(beta21.gva[6, , index0])),
    rep("5*10^4", length(beta21.gva[7, , index0])),
    rep("9*10^4", length(beta21.gva[8, , index0]))
  ),
  value = c(
    beta21.gva[1, , index0],
    beta21.gva[2, , index0],
    beta21.gva[3, , index0],
    beta21.gva[4, , index0],
    beta21.gva[5, , index0],
    beta21.gva[6, , index0],
    beta21.gva[7, , index0],
    beta21.gva[8, , index0]
  )
)
data_gva_beta21$group <- factor(data_gva_beta21$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta21_gva <- ggplot(data_gva_beta21, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[21])), title = expression("GVA: Boxplots of " * hat(beta)[21]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.25, -0.75) 
data_gva_beta22 <- data.frame(
  group = c(
    rep("1*10^2", length(beta22.gva[1, , index0])),
    rep("5*10^2", length(beta22.gva[2, , index0])),
    rep("1*10^3", length(beta22.gva[3, , index0])),
    rep("5*10^3", length(beta22.gva[4, , index0])),
    rep("1*10^4", length(beta22.gva[5, , index0])),
    rep("3*10^4", length(beta22.gva[6, , index0])),
    rep("5*10^4", length(beta22.gva[7, , index0])),
    rep("9*10^4", length(beta22.gva[8, , index0]))
  ),
  value = c(
    beta22.gva[1, , index0],
    beta22.gva[2, , index0],
    beta22.gva[3, , index0],
    beta22.gva[4, , index0],
    beta22.gva[5, , index0],
    beta22.gva[6, , index0],
    beta22.gva[7, , index0],
    beta22.gva[8, , index0]
  )
)
data_gva_beta22$group <- factor(data_gva_beta22$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta22_gva <- ggplot(data_gva_beta22, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[22])), title = expression("GVA: Boxplots of " * hat(beta)[22]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.25, -0.75) 
data_gva_sigma1 <- data.frame(
  group = c(
    rep("1*10^2", length(sigma1.gva[1, , index0])),
    rep("5*10^2", length(sigma1.gva[2, , index0])),
    rep("1*10^3", length(sigma1.gva[3, , index0])),
    rep("5*10^3", length(sigma1.gva[4, , index0])),
    rep("1*10^4", length(sigma1.gva[5, , index0])),
    rep("3*10^4", length(sigma1.gva[6, , index0])),
    rep("5*10^4", length(sigma1.gva[7, , index0])),
    rep("9*10^4", length(sigma1.gva[8, , index0]))
  ),
  value = c(
    sigma1.gva[1, , index0],
    sigma1.gva[2, , index0],
    sigma1.gva[3, , index0],
    sigma1.gva[4, , index0],
    sigma1.gva[5, , index0],
    sigma1.gva[6, , index0],
    sigma1.gva[7, , index0],
    sigma1.gva[8, , index0]
  )
)
data_gva_sigma1$group <- factor(data_gva_sigma1$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_sigma1_gva <- ggplot(data_gva_sigma1, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(sigma[1])), title = expression("GVA: Boxplots of " * hat(sigma)[1]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(0, 2) 
data_gva_sigma2 <- data.frame(
  group = c(
    rep("1*10^2", length(sigma2.gva[1, , index0])),
    rep("5*10^2", length(sigma2.gva[2, , index0])),
    rep("1*10^3", length(sigma2.gva[3, , index0])),
    rep("5*10^3", length(sigma2.gva[4, , index0])),
    rep("1*10^4", length(sigma2.gva[5, , index0])),
    rep("3*10^4", length(sigma2.gva[6, , index0])),
    rep("5*10^4", length(sigma2.gva[7, , index0])),
    rep("9*10^4", length(sigma2.gva[8, , index0]))
  ),
  value = c(
    sigma2.gva[1, , index0],
    sigma2.gva[2, , index0],
    sigma2.gva[3, , index0],
    sigma2.gva[4, , index0],
    sigma2.gva[5, , index0],
    sigma2.gva[6, , index0],
    sigma2.gva[7, , index0],
    sigma2.gva[8, , index0]
  )
)
data_gva_sigma2$group <- factor(data_gva_sigma2$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_sigma2_gva <- ggplot(data_gva_sigma2, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size", y = expression(hat(sigma[2])), title = expression("GVA: Boxplots of " * hat(sigma)[2]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(0, 2) 
data_gva_rho <- data.frame(
  group = c(
    rep("1*10^2", length(rho.gva[1, , index0])),
    rep("5*10^2", length(rho.gva[2, , index0])),
    rep("1*10^3", length(rho.gva[3, , index0])),
    rep("5*10^3", length(rho.gva[4, , index0])),
    rep("1*10^4", length(rho.gva[5, , index0])),
    rep("3*10^4", length(rho.gva[6, , index0])),
    rep("5*10^4", length(rho.gva[7, , index0])),
    rep("9*10^4", length(rho.gva[8, , index0]))
  ),
  value = c(
    rho.gva[1, , index0],
    rho.gva[2, , index0],
    rho.gva[3, , index0],
    rho.gva[4, , index0],
    rho.gva[5, , index0],
    rho.gva[6, , index0],
    rho.gva[7, , index0],
    rho.gva[8, , index0]
  )
)
data_gva_rho$group <- factor(data_gva_rho$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_rho_gva <- ggplot(data_gva_rho, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = 0.56, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(rho)), title = expression("GVA: Boxplots of " * hat(rho)))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1, 1.4) 
#mcmc_boxplot
data_mcmc_beta01 <- data.frame(
  group = c(
    rep("1*10^2", length(beta01.mcmc[1, , index0])),
    rep("5*10^2", length(beta01.mcmc[2, , index0])),
    rep("1*10^3", length(beta01.mcmc[3, , index0])),
    rep("5*10^3", length(beta01.mcmc[4, , index0])),
    rep("1*10^4", length(beta01.mcmc[5, , index0])),
    rep("3*10^4", length(beta01.mcmc[6, , index0])),
    rep("5*10^4", length(beta01.mcmc[7, , index0])),
    rep("9*10^4", length(beta01.mcmc[8, , index0]))
  ),
  value = c(
    beta01.mcmc[1, , index0],
    beta01.mcmc[2, , index0],
    beta01.mcmc[3, , index0],
    beta01.mcmc[4, , index0],
    beta01.mcmc[5, , index0],
    beta01.mcmc[6, , index0],
    beta01.mcmc[7, , index0],
    beta01.mcmc[8, , index0]
  )
)
data_mcmc_beta01$group <- factor(data_mcmc_beta01$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta01_mcmc <- ggplot(data_mcmc_beta01, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[0*1])), title = expression("MCMC: Boxplots of " * hat(beta)[0*1]))+
  theme_minimal() +
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-2, 0) 
data_mcmc_beta11 <- data.frame(
  group = c(
    rep("1*10^2", length(beta11.mcmc[1, , index0])),
    rep("5*10^2", length(beta11.mcmc[2, , index0])),
    rep("1*10^3", length(beta11.mcmc[3, , index0])),
    rep("5*10^3", length(beta11.mcmc[4, , index0])),
    rep("1*10^4", length(beta11.mcmc[5, , index0])),
    rep("3*10^4", length(beta11.mcmc[6, , index0])),
    rep("5*10^4", length(beta11.mcmc[7, , index0])),
    rep("9*10^4", length(beta11.mcmc[8, , index0]))
  ),
  value = c(
    beta11.mcmc[1, , index0],
    beta11.mcmc[2, , index0],
    beta11.mcmc[3, , index0],
    beta11.mcmc[4, , index0],
    beta11.mcmc[5, , index0],
    beta11.mcmc[6, , index0],
    beta11.mcmc[7, , index0],
    beta11.mcmc[8, , index0]
  )
)
data_mcmc_beta11$group <- factor(data_mcmc_beta11$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta11_mcmc <- ggplot(data_mcmc_beta11, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[11])), title = expression("MCMC: Boxplots of " * hat(beta)[11]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.5, -0.5) 
data_mcmc_beta12 <- data.frame(
  group = c(
    rep("1*10^2", length(beta12.mcmc[1, , index0])),
    rep("5*10^2", length(beta12.mcmc[2, , index0])),
    rep("1*10^3", length(beta12.mcmc[3, , index0])),
    rep("5*10^3", length(beta12.mcmc[4, , index0])),
    rep("1*10^4", length(beta12.mcmc[5, , index0])),
    rep("3*10^4", length(beta12.mcmc[6, , index0])),
    rep("5*10^4", length(beta12.mcmc[7, , index0])),
    rep("9*10^4", length(beta12.mcmc[8, , index0]))
  ),
  value = c(
    beta12.mcmc[1, , index0],
    beta12.mcmc[2, , index0],
    beta12.mcmc[3, , index0],
    beta12.mcmc[4, , index0],
    beta12.mcmc[5, , index0],
    beta12.mcmc[6, , index0],
    beta12.mcmc[7, , index0],
    beta12.mcmc[8, , index0]
  )
)
data_mcmc_beta12$group <- factor(data_mcmc_beta12$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta12_mcmc <- ggplot(data_mcmc_beta12, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[12])), title = expression("MCMC: Boxplots of " * hat(beta)[12]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.25, -0.75) 
data_mcmc_beta02 <- data.frame(
  group = c(
    rep("1*10^2", length(beta02.mcmc[1, , index0])),
    rep("5*10^2", length(beta02.mcmc[2, , index0])),
    rep("1*10^3", length(beta02.mcmc[3, , index0])),
    rep("5*10^3", length(beta02.mcmc[4, , index0])),
    rep("1*10^4", length(beta02.mcmc[5, , index0])),
    rep("3*10^4", length(beta02.mcmc[6, , index0])),
    rep("5*10^4", length(beta02.mcmc[7, , index0])),
    rep("9*10^4", length(beta02.mcmc[8, , index0]))
  ),
  value = c(
    beta02.mcmc[1, , index0],
    beta02.mcmc[2, , index0],
    beta02.mcmc[3, , index0],
    beta02.mcmc[4, , index0],
    beta02.mcmc[5, , index0],
    beta02.mcmc[6, , index0],
    beta02.mcmc[7, , index0],
    beta02.mcmc[8, , index0]
  )
)
data_mcmc_beta02$group <- factor(data_mcmc_beta02$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta02_mcmc <- ggplot(data_mcmc_beta02, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[0*2])), title = expression("MCMC: Boxplots of " * hat(beta)[0*2]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-2, 0) 
data_mcmc_beta21 <- data.frame(
  group = c(
    rep("1*10^2", length(beta21.mcmc[1, , index0])),
    rep("5*10^2", length(beta21.mcmc[2, , index0])),
    rep("1*10^3", length(beta21.mcmc[3, , index0])),
    rep("5*10^3", length(beta21.mcmc[4, , index0])),
    rep("1*10^4", length(beta21.mcmc[5, , index0])),
    rep("3*10^4", length(beta21.mcmc[6, , index0])),
    rep("5*10^4", length(beta21.mcmc[7, , index0])),
    rep("9*10^4", length(beta21.mcmc[8, , index0]))
  ),
  value = c(
    beta21.mcmc[1, , index0],
    beta21.mcmc[2, , index0],
    beta21.mcmc[3, , index0],
    beta21.mcmc[4, , index0],
    beta21.mcmc[5, , index0],
    beta21.mcmc[6, , index0],
    beta21.mcmc[7, , index0],
    beta21.mcmc[8, , index0]
  )
)
data_mcmc_beta21$group <- factor(data_mcmc_beta21$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta21_mcmc <- ggplot(data_mcmc_beta21, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[21])), title = expression("MCMC: Boxplots of " * hat(beta)[21]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.25, -0.75) 
data_mcmc_beta22 <- data.frame(
  group = c(
    rep("1*10^2", length(beta22.mcmc[1, , index0])),
    rep("5*10^2", length(beta22.mcmc[2, , index0])),
    rep("1*10^3", length(beta22.mcmc[3, , index0])),
    rep("5*10^3", length(beta22.mcmc[4, , index0])),
    rep("1*10^4", length(beta22.mcmc[5, , index0])),
    rep("3*10^4", length(beta22.mcmc[6, , index0])),
    rep("5*10^4", length(beta22.mcmc[7, , index0])),
    rep("9*10^4", length(beta22.mcmc[8, , index0]))
  ),
  value = c(
    beta22.mcmc[1, , index0],
    beta22.mcmc[2, , index0],
    beta22.mcmc[3, , index0],
    beta22.mcmc[4, , index0],
    beta22.mcmc[5, , index0],
    beta22.mcmc[6, , index0],
    beta22.mcmc[7, , index0],
    beta22.mcmc[8, , index0]
  )
)
data_mcmc_beta22$group <- factor(data_mcmc_beta22$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_beta22_mcmc <- ggplot(data_mcmc_beta22, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = -1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(beta[22])), title = expression("MCMC: Boxplots of " * hat(beta)[22]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1.25, -0.75) 
data_mcmc_sigma1 <- data.frame(
  group = c(
    rep("1*10^2", length(sigma1.mcmc[1, , index0])),
    rep("5*10^2", length(sigma1.mcmc[2, , index0])),
    rep("1*10^3", length(sigma1.mcmc[3, , index0])),
    rep("5*10^3", length(sigma1.mcmc[4, , index0])),
    rep("1*10^4", length(sigma1.mcmc[5, , index0])),
    rep("3*10^4", length(sigma1.mcmc[6, , index0])),
    rep("5*10^4", length(sigma1.mcmc[7, , index0])),
    rep("9*10^4", length(sigma1.mcmc[8, , index0]))
  ),
  value = c(
    sigma1.mcmc[1, , index0],
    sigma1.mcmc[2, , index0],
    sigma1.mcmc[3, , index0],
    sigma1.mcmc[4, , index0],
    sigma1.mcmc[5, , index0],
    sigma1.mcmc[6, , index0],
    sigma1.mcmc[7, , index0],
    sigma1.mcmc[8, , index0]
  )
)
data_mcmc_sigma1$group <- factor(data_mcmc_sigma1$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_sigma1_mcmc <- ggplot(data_mcmc_sigma1, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(sigma[1])), title = expression("MCMC: Boxplots of " * hat(sigma)[1]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(0, 2) 
data_mcmc_sigma2 <- data.frame(
  group = c(
    rep("1*10^2", length(sigma2.mcmc[1, , index0])),
    rep("5*10^2", length(sigma2.mcmc[2, , index0])),
    rep("1*10^3", length(sigma2.mcmc[3, , index0])),
    rep("5*10^3", length(sigma2.mcmc[4, , index0])),
    rep("1*10^4", length(sigma2.mcmc[5, , index0])),
    rep("3*10^4", length(sigma2.mcmc[6, , index0])),
    rep("5*10^4", length(sigma2.mcmc[7, , index0])),
    rep("9*10^4", length(sigma2.mcmc[8, , index0]))
  ),
  value = c(
    sigma2.mcmc[1, , index0],
    sigma2.mcmc[2, , index0],
    sigma2.mcmc[3, , index0],
    sigma2.mcmc[4, , index0],
    sigma2.mcmc[5, , index0],
    sigma2.mcmc[6, , index0],
    sigma2.mcmc[7, , index0],
    sigma2.mcmc[8, , index0]
  )
)
data_mcmc_sigma2$group <- factor(data_mcmc_sigma2$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_sigma2_mcmc <- ggplot(data_mcmc_sigma2, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(sigma[2])), title = expression("MCMC: Boxplots of " * hat(sigma)[2]))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(0, 2) 
data_mcmc_rho <- data.frame(
  group = c(
    rep("1*10^2", length(rho.mcmc[1, , index0])),
    rep("5*10^2", length(rho.mcmc[2, , index0])),
    rep("1*10^3", length(rho.mcmc[3, , index0])),
    rep("5*10^3", length(rho.mcmc[4, , index0])),
    rep("1*10^4", length(rho.mcmc[5, , index0])),
    rep("3*10^4", length(rho.mcmc[6, , index0])),
    rep("5*10^4", length(rho.mcmc[7, , index0])),
    rep("9*10^4", length(rho.mcmc[8, , index0]))
  ),
  value = c(
    rho.mcmc[1, , index0],
    rho.mcmc[2, , index0],
    rho.mcmc[3, , index0],
    rho.mcmc[4, , index0],
    rho.mcmc[5, , index0],
    rho.mcmc[6, , index0],
    rho.mcmc[7, , index0],
    rho.mcmc[8, , index0]
  )
)
data_mcmc_rho$group <- factor(data_mcmc_rho$group,levels = c("1*10^2","5*10^2","1*10^3","5*10^3","1*10^4","3*10^4","5*10^4","9*10^4"))
boxplot_rho_mcmc <- ggplot(data_mcmc_rho, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = "gray",outlier.color = "gray",outlier.fill = "gray") +
  geom_hline(yintercept = 0.56, color = "red", linetype = "solid", size = 0.8)+
  labs(x = "Sample Size ", y = expression(hat(rho)), title = expression("MCMC: Boxplots of " * hat(rho)))+
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)  # 居中标题
  ) +
  ylim(-1,1.4) 
# 保存当前的 R 会话
save.image("C:/Users/86188/Desktop/simulation for multi Poisson response/my_workspace.RData")
