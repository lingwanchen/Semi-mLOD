nega<-function(x) -x
inv.nega<-function(x) -x
exp.ne<-function(x) exp(-x)
inv.exp.ne<-function(x) -log(x)

cummatrix = function(ee){ 
  n2 <- dim(ee)[1]
  n4 <- dim(ee)[2]
  
  for (i in 1:n4) ee[,i] <- cumsum(ee[,i])
  for (i in 1:n2) ee[i,] <- cumsum(ee[i,])
  
  return(ee)
}

cumprodmatrix = function(ee){
  n2 <- dim(ee)[1]
  n4 <- dim(ee)[2]
  
  for (i in 1:n4) ee[,i] <- cumprod(ee[,i])
  for (i in 1:n2) ee[i,] <- cumprod(ee[i,])
  
  return(ee)
}

survf<-function(time,delta){
  ndim <- length(time)        
  uniquetime1 <- unique(sort(time))
  npts1 <- length(uniquetime1)
  
  tmp1 <- matrix(rep(time,each = npts1),npts1,ndim) 
  tmp3 <- matrix(rep(delta,each = npts1),npts1,ndim) 
  
  lambda1 <- apply(ifelse(tmp1==uniquetime1 & tmp3 ==1,1,0),1,sum)/apply(ifelse(tmp1>=uniquetime1,1,0),1,sum)   
  S <- cumprod(1-lambda1)
  jump <- -diff(c(1,S))
  return(list(time=uniquetime1, S=S, prob=jump))
}

DB.BS<-function(DD){
  ndim <- dim(DD)[1]        
  uniquetime1 <- unique(sort(DD[,1]))  # for Y
  uniquetime2 <- unique(sort(DD[,2]))  # for T
  
  npts1 <- length(uniquetime1)
  npts2 <- length(uniquetime2)
  
  tmp1 <- matrix(rep(DD[,1],each = npts1),npts1,ndim) # for Y 
  tmp3 <- matrix(rep(DD[,3],each = npts1),npts1,ndim) # for delta_Y
  
  tmp2 <- matrix(rep(DD[,2],each = npts2),npts2,ndim) # for T
  tmp4 <- matrix(rep(DD[,4],each = npts2),npts2,ndim) # for delta_T
  
  numrisk <- ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))
  lambda11 <- ifelse(tmp1==uniquetime1 & tmp3!=0,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 !=0,1,0))/numrisk
  
  lambda10.left <- ifelse(tmp1==uniquetime1 & tmp3 !=0,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))/numrisk
  lambda01.up <- ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 !=0,1,0))/numrisk
  
  numrisk <- NULL
  
  lambda10.left <- ifelse(lambda10.left == "NaN", 0, lambda10.left)
  lambda01.up <- ifelse(lambda01.up == "NaN", 0, lambda01.up)
  lambda11 <- ifelse(lambda11 == "NaN", 0, lambda11)   
  
  # for Y
  lambda1 <- apply(ifelse(tmp1==uniquetime1 & tmp3 != 0,1,0),1,sum)/apply(ifelse(tmp1>=uniquetime1,1,0),1,sum)   
  S1 <- cumprod(1-lambda1)
  
  # for T
  lambda2 <- apply(ifelse(tmp2==uniquetime2 & tmp4 != 0,1,0),1,sum)/apply(ifelse(tmp2>=uniquetime2,1,0),1,sum)   
  S2 <- cumprod(1-lambda2)
  
  ## Dabrowska estimator for the bivariate survival function
  id11.00 <- apply(ifelse(tmp1==uniquetime1 & tmp3==0,1,0),1,sum)
  id21.00 <- apply(ifelse(tmp2==uniquetime2 & tmp4==0,1,0),1,sum)
  
  tmp1 = tmp2 = tmp3 = tmp4 = NULL
  
  # factor
  
  factor <- ifelse(lambda10.left!=1 & lambda01.up!=1, 1-(lambda10.left*lambda01.up-lambda11)/(1-lambda10.left)/(1-lambda01.up), 1)
  factor <- cumprodmatrix(factor)   
  factor <- ifelse(is.na(factor),0,factor)
  
  SS <- S1%*%t(S2)*factor
  
  factor <- NULL
  
  # to handle non-monotonicity of Dabrowska's estimator, use cummin
  SS <- apply(SS,2,cummin)
  SS <- t(apply(SS,1,cummin))
  
  return(list(BS=SS[(id11.00==0),(id21.00==0)],
              S1=S1[(id11.00==0)],
              S2=S2[(id21.00==0)],
              t1=uniquetime1[(id11.00==0)],
              t2=uniquetime2[(id21.00==0)]))
}

BSF2<-function(DD){
  ndim <- dim(DD)[1]        
  uniquetime1 <- unique(sort(DD[(DD[,3]==1),1]))  # for Y
  uniquetime2 <- unique(sort(DD[(DD[,4]==1),2]))  # for T
  
  npts1 <- length(uniquetime1)
  npts2 <- length(uniquetime2)
  
  tmp1 <- matrix(rep(DD[,1],each = npts1),npts1,ndim) # for Y 
  tmp3 <- matrix(rep(DD[,3],each = npts1),npts1,ndim) # for delta_Y
  
  tmp2 <- matrix(rep(DD[,2],each = npts2),npts2,ndim) # for T
  tmp4 <- matrix(rep(DD[,4],each = npts2),npts2,ndim) # for delta_T
  
  numrisk <- ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))
  lambda11 <- ifelse(tmp1==uniquetime1 & tmp3!=0,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 !=0,1,0))/numrisk
  
  lambda10.left <- ifelse(tmp1==uniquetime1 & tmp3 !=0,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))/numrisk
  lambda01.up <- ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 !=0,1,0))/numrisk
  
  numrisk <- NULL
  
  lambda10.left <- ifelse(lambda10.left == "NaN", 0, lambda10.left)
  lambda01.up <- ifelse(lambda01.up == "NaN", 0, lambda01.up)
  lambda11 <- ifelse(lambda11 == "NaN", 0, lambda11)   
  
  # for Y
  lambda1 <- apply(ifelse(tmp1==uniquetime1 & tmp3 != 0,1,0),1,sum)/apply(ifelse(tmp1>=uniquetime1,1,0),1,sum)   
  S1 <- cumprod(1-lambda1)
  
  # for T
  lambda2 <- apply(ifelse(tmp2==uniquetime2 & tmp4 != 0,1,0),1,sum)/apply(ifelse(tmp2>=uniquetime2,1,0),1,sum)   
  S2 <- cumprod(1-lambda2)
  
  ## Dabrowska estimator for the bivariate survival function
  
  tmp1 = tmp2 = tmp3 = tmp4 = NULL
  
  # factor
  
  factor <- ifelse(lambda10.left!=1 & lambda01.up!=1, 1-(lambda10.left*lambda01.up-lambda11)/(1-lambda10.left)/(1-lambda01.up), 1)
  factor <- cumprodmatrix(factor)   
  factor <- ifelse(is.na(factor),0,factor)
  
  SS <- S1%*%t(S2)*factor
  
  factor <- NULL
  
  # to handle non-monotonicity of Dabrowska's estimator, use cummin
  SS <- apply(SS,2,cummin)
  SS <- t(apply(SS,1,cummin))
  
  return(list(BS=SS,S1=S1,S2=S2,t1=uniquetime1,t2=uniquetime2 ))
}

TriKM <- function(DD){
  n <- dim(DD)[1]
  I <- length(unique(DD[(DD[,4]==1),1])) #sum(DD[,4])
  J <- length(unique(DD[(DD[,5]==1),2])) #sum(DD[,5])
  K <- length(unique(DD[(DD[,6]==1),3])) #sum(DD[,6])
  
  Fhat <- array(rep(0, (I+1)*(J+1)*(K+1)), dim=c(I+1, J+1, K+1))
  Fhat[1, 1, 1] <- 1
  
  ##### Marginal & Bivariate #####
  D12<-DD[,c(1,2,4,5)]
  D13<-DD[,c(1,3,4,6)]
  D23<-DD[,c(2,3,5,6)]
  
  BS12<-BSF2(D12)
  BS13<-BSF2(D13)
  BS23<-BSF2(D23)
  
  S1 <- BS12$S1
  S2 <- BS12$S2
  S3 <- BS13$S2
  
  S12 <- BS12$BS
  S13 <- BS13$BS
  S23 <- BS23$BS
  
  BS12 <- BS13 <- BS23 <- NULL
  D12 <- D23 <- D13 <-NULL
  
  Fhat[-1, 1, 1] <- S1
  Fhat[1, -1, 1] <- S2
  Fhat[1, 1, -1] <- S3
  Fhat[-1, -1, 1] <- S12
  Fhat[-1, 1, -1] <- S13
  Fhat[1, -1, -1] <- S23
  
  ##### Trivariate #####
  Y1 <- DD[,1]
  Y2 <- DD[,2]
  Y3 <- DD[,3]
  Delta1 <- DD[,4]
  Delta2 <- DD[,5]
  Delta3 <- DD[,6]
  
  T1 <- unique(sort(Y1[Delta1==1]))
  T2 <- unique(sort(Y2[Delta2==1]))
  T3 <- unique(sort(Y3[Delta3==1]))
  
  a0<-Sys.time()
  for (i in 1:I){
    
    for (j in 1:J){
      
      for (k in 1:K){
        
        if (Fhat[i+1, j+1, k]*Fhat[i+1, j, k+1]*Fhat[i, j+1, k+1]==0){
          
          Fhat[i+1, j+1, k+1] <- 0
          
        }else{
          
          rijk <- sum(Y1>=T1[i] & Y2>=T2[j] & Y3>=T3[k])
          
          if (rijk>0){
            
            Ftilda_111 <- sum(Y1>T1[i] & Y2>T2[j] & Y3>T3[k])/rijk
            Ftilda_110 <- sum(Y1>T1[i] & Y2>T2[j] & Y3>=T3[k])/rijk
            Ftilda_101 <- sum(Y1>T1[i] & Y2>=T2[j] & Y3>T3[k])/rijk
            Ftilda_011 <- sum(Y1>=T1[i] & Y2>T2[j] & Y3>T3[k])/rijk
            Ftilda_100 <- sum(Y1>T1[i] & Y2>=T2[j] & Y3>=T3[k])/rijk
            Ftilda_010 <- sum(Y1>=T1[i] & Y2>T2[j] & Y3>=T3[k])/rijk
            Ftilda_001 <- sum(Y1>=T1[i] & Y2>=T2[j] & Y3>T3[k])/rijk
            
            if (Ftilda_110*Ftilda_101*Ftilda_011*Ftilda_001*Ftilda_010*Ftilda_100>0){
              
              Fhat[i+1, j+1, k+1] <- Fhat[i+1, j+1, k]*Fhat[i+1, j, k+1]*Fhat[i, j+1, k+1]*Fhat[i, j, k]/
                (Fhat[i+1, j, k]*Fhat[i, j+1, k]*Fhat[i, j, k+1])*Ftilda_111*Ftilda_100*Ftilda_010*Ftilda_001/(Ftilda_110*Ftilda_101*Ftilda_011)
            }else{
              Fhat[i+1, j+1, k+1] <- Fhat[i+1, j+1, k]*Fhat[i+1, j, k+1]*Fhat[i, j+1, k+1]*Fhat[i, j, k]/(Fhat[i+1, j, k]*Fhat[i, j+1, k]*Fhat[i, j, k+1])
            }
          }
          else{
            Fhat[i+1, j+1, k+1] <- Fhat[i+1, j+1, k]*Fhat[i+1, j, k+1]*Fhat[i, j+1, k+1]*Fhat[i, j, k]/(Fhat[i+1, j, k]*Fhat[i, j+1, k]*Fhat[i, j, k+1])
          }
          
        }
        
      }
      
    }
    
  }
  a1<-Sys.time()
  
  return(list(TSF=Fhat, S1=Fhat[, 1, 1], S2=Fhat[1, , 1], S3=Fhat[1, 1, ],
              T1=T1, T2=T2, T3=T3,
              S12=Fhat[, , 1], S13=Fhat[, 1, ], S23=Fhat[1, , ]))
}

coef.rankAFT=function (x, y, delta){
  ynew <- 1000 * (length(y))^2
  data1 <- data.frame(y, x)
  options(contrasts = c("contr.treatment", "contr.poly"))
  tempfit <- lm(y ~ ., x = TRUE, y = TRUE, data = data1)
  x <- as.matrix(tempfit$x[, -1])
  xn <- dimnames(x)[[2]]
  y <- tempfit$y
  dimnum <- dim(x)
  n1 <- dimnum[1]
  n2 <- dimnum[2]
  yy0 <- rep(y, rep(n1, n1))
  delta1 <- rep(delta, rep(n1, n1))
  yy1 <- rep(y, n1)
  yy <- delta1 * (yy0 - yy1)
  xx0 <- matrix(rep(as.vector(x), rep(n1, n1 * n2)), nrow = n1 * n1)
  xx1 <- t(matrix(rep(as.vector(t(x)), n1), nrow = n2))
  xx <- xx0 - xx1
  xxdif <- xx * delta1
  xnew <- apply(xxdif, 2, sum)
  xnew <- rbind(xxdif, -xnew)
  yynew <- c(yy, ynew)
  dimnames(xnew) <- list(NULL, xn)
  beta <- rq(yynew ~ xnew - 1, method = "fn")$coef
  #beta <- fit$coef
  
  res <- y-as.matrix(x)%*%beta
  # K-M estimate of the survival function of the error term (at the estimated beta)
  KM.fit <- survf(res,delta)
  # KM.fit$time is ordered
  time <- unique(sort(res))
  #surv <- KM.fit$sur # estimated survival function for the residual error
  jump <- KM.fit$prob #-surv # drop the last one if it is not a failure (in that case, surv[q-1]=surv[q])
  #jump <- c(1,surv[-q])-c(surv[-q],0) #force last one to be failure if it is not a failure
  alpha <- time%*%jump # intercept estimator of the AFT model
  
  #  predmatrix <- x - t(matrix(rep(apply(x, 2, mean), n1), ncol = n1))
  #  residuals <- y - predmatrix %*% as.matrix(betag)
  
  return(c(alpha, beta))
}


pslogL.LM.p1<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  sigma.0<-parms$sigma.0
  
  err.y<-data$Y-theta%*%t(XZ)
  
  if(length(idx.1)>0L){
    pl1.sum<-err.y[idx.1]%*%XZ[idx.1,]/n
  } else{ pl1.sum <- 0 }
  
  pl2.0<-matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.0]==n1)
  if( length(rm.i)>0) idx.0 <- idx.0[-rm.i]
  
  if(length(idx.0)>0L){
    for (j in idx.0){
      idx.h <- ord.res[j]
      prob<-ff[(idx.h+1):n1]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,])
      XZ1 <- cbind(x.s, z1.hat)
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl2.0<-rbind(pl2.0, nfyz1/dfyz1)
      }
    }
    pl2.sum<-colSums(pl2.0)/n
  } else{ pl2.sum <- 0 }
  
  pl1.sum+pl2.sum
}

pslogL.BN.p1<-function(theta, parms){
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(1/(1+exp(-ymu)))
  
  if(length(idx.1)>0L){
    pl1.sum<-err.y[idx.1]%*%XZ[idx.1,]/n
  } else{ pl1.sum <- 0 }
  
  pl2.0<-matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.0]==n1)
  if( length(rm.i)>0) idx.0 <- idx.0[-rm.i]
  
  if(length(idx.0)>0L){
    for (j in idx.0){
      idx.h <- ord.res[j]
      prob<-ff[(idx.h+1):n1]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,])
      XZ1 <- cbind(x.s, z1.hat)
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl2.0<-rbind(pl2.0, nfyz1/dfyz1)
      }
    }
    pl2.sum<-colSums(pl2.0)/n
  } else{ pl2.sum <- 0 }
  
  pl1.sum+pl2.sum
}

pslogL.PO.p1<-function(theta, parms){
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(exp(ymu))
  
  if(length(idx.1)>0L){
    pl1.sum<-err.y[idx.1]%*%XZ[idx.1,]/n
  } else{ pl1.sum <- 0 }
  
  pl2.0<-matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.0]==n1)
  if( length(rm.i)>0) idx.0 <- idx.0[-rm.i]
  
  if(length(idx.0)>0L){
    for (j in idx.0){
      idx.h <- ord.res[j]
      prob<-ff[(idx.h+1):n1]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,])
      XZ1 <- cbind(x.s, z1.hat)
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl2.0<-rbind(pl2.0, nfyz1/dfyz1)
      }
    }
    pl2.sum<-colSums(pl2.0)/n
  } else{ pl2.sum <- 0 }
  
  pl1.sum+pl2.sum
}

pslogL.LM.p2<-function(theta, parms){
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  sigma.0<-parms$sigma.0
  
  err.y<-data$Y-theta%*%t(XZ)

  if(length(idx.11)>0L){
    pl1.sum<-err.y[idx.11]%*%XZ[idx.11,]/n
  } else{ pl1.sum <- 0 }
  
  
  pl2.0 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.01,1]==n1)
  if( length(rm.i)>0) idx.01 <- idx.01[-rm.i]
  
  if(length(idx.01)>0L){
    for (j in idx.01){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      
      prob<-ff[(idx.h+1):n1, (idx.v)]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j])
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl2.0<-rbind(pl2.0, nfyz1/dfyz1)
      }
    }
    pl2.sum<-colSums(pl2.0)/n
    
  } else{ pl2.sum <- 0 }
  
  pl3.0 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.10,2]==n2)
  if( length(rm.i)>0) idx.10 <- idx.10[-rm.i]
  
  if(length(idx.10)>0L){
    for (j in idx.10){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      
      prob<-ff[(idx.h), (idx.v+1):n2]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,])
      XZ2 <- cbind(x.s, data$z1st[j], z2.hat)
      
      err.yz2 <- c(data$Y[j] - theta%*%t(XZ2))
      fyz2 <- c(dnorm(err.yz2, sd=sigma.0))
      nfyz2 <- c(fyz2*prob*err.yz2)%*%XZ2
      dfyz2 <- sum(fyz2*prob)
      if(dfyz2>0L){
        pl3.0<-rbind(pl3.0, nfyz2/dfyz2)
      }
    }
    
    pl3.sum<-colSums(pl3.0)/n
  } else{ pl3.sum <- 0 }
  
  pl4.01 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.00,1]==n1)
  rm.2i <- which(ord.res[idx.00,2]==n2)
  
  if( (length(rm.1i)+length(rm.2i))>0) idx.00 <- idx.00[-unique(c(rm.1i,rm.2i))]
  
  if(length(idx.00)>0L){
    for (j in idx.00){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      prob<-ff[(idx.h+1):n1, (idx.v+1):n2]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1))
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl4.01<-rbind(pl4.01, nfyz12/dfyz12)
      }
    }
    
    pl4.sum<-colSums(pl4.01)/n
  } else{ pl4.sum <- 0 }
  
  pl1.sum+pl2.sum+pl3.sum+pl4.sum
  
}

pslogL.BN.p2<-function(theta, parms){
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(1/(1+exp(-ymu)))
  
  if(length(idx.11)>0L){
    pl1.sum<-err.y[idx.11]%*%XZ[idx.11,]/n
  } else{ pl1.sum <- 0 }
  
  
  pl2.0 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.01,1]==n1)
  if( length(rm.i)>0) idx.01 <- idx.01[-rm.i]
  
  if(length(idx.01)>0L){
    for (j in idx.01){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      
      prob<-ff[(idx.h+1):n1, (idx.v)]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-data$Y[j]- y.pr
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl2.0<-rbind(pl2.0, nfyz1/dfyz1)
      }
    }
    pl2.sum<-colSums(pl2.0)/n
    
  } else{ pl2.sum <- 0 }
  
  pl3.0 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.10,2]==n2)
  if( length(rm.i)>0) idx.10 <- idx.10[-rm.i]
  
  if(length(idx.10)>0L){
    for (j in idx.10){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      
      prob<-ff[(idx.h), (idx.v+1):n2]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,])
      XZ2 <- cbind(x.s, data$z1st[j], z2.hat)
      
      y.mu <- theta%*%t(XZ2)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz2<-data$Y[j]- y.pr
      
      fyz2 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz2 <- c(fyz2*prob*err.yz2)%*%XZ2
      dfyz2 <- sum(fyz2*prob)
      if(dfyz2>0L){
        pl3.0<-rbind(pl3.0, nfyz2/dfyz2)
      }
    }
    
    pl3.sum<-colSums(pl3.0)/n
  } else{ pl3.sum <- 0 }
  
  pl4.01 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.00,1]==n1)
  rm.2i <- which(ord.res[idx.00,2]==n2)
  
  if( (length(rm.1i)+length(rm.2i))>0) idx.00 <- idx.00[-unique(c(rm.1i,rm.2i))]
  
  if(length(idx.00)>0L){
    for (j in idx.00){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      prob<-ff[(idx.h+1):n1, (idx.v+1):n2]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl4.01<-rbind(pl4.01, nfyz12/dfyz12)
      }
    }
    
    pl4.sum<-colSums(pl4.01)/n
  } else{ pl4.sum <- 0 }
  
  pl1.sum+pl2.sum+pl3.sum+pl4.sum
  
}

pslogL.PO.p2<-function(theta, parms){
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(exp(ymu))
  
  if(length(idx.11)>0L){
    pl1.sum<-err.y[idx.11]%*%XZ[idx.11,]/n
  } else{ pl1.sum <- 0 }
  
  
  pl2.0 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.01,1]==n1)
  if( length(rm.i)>0) idx.01 <- idx.01[-rm.i]
  
  if(length(idx.01)>0L){
    for (j in idx.01){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      
      prob<-ff[(idx.h+1):n1, (idx.v)]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-data$Y[j]- y.pr
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl2.0<-rbind(pl2.0, nfyz1/dfyz1)
      }
    }
    pl2.sum<-colSums(pl2.0)/n
    
  } else{ pl2.sum <- 0 }
  
  pl3.0 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.10,2]==n2)
  if( length(rm.i)>0) idx.10 <- idx.10[-rm.i]
  
  if(length(idx.10)>0L){
    for (j in idx.10){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      
      prob<-ff[(idx.h), (idx.v+1):n2]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,])
      XZ2 <- cbind(x.s, data$z1st[j], z2.hat)
      
      y.mu <- theta%*%t(XZ2)
      y.pr <- c(exp(y.mu))
      err.yz2<-data$Y[j]- y.pr
      fyz2 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz2 <- c(fyz2*prob*err.yz2)%*%XZ2
      dfyz2 <- sum(fyz2*prob)
      if(dfyz2>0L){
        pl3.0<-rbind(pl3.0, nfyz2/dfyz2)
      }
    }
    
    pl3.sum<-colSums(pl3.0)/n
  } else{ pl3.sum <- 0 }
  
  pl4.01 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.00,1]==n1)
  rm.2i <- which(ord.res[idx.00,2]==n2)
  
  if( (length(rm.1i)+length(rm.2i))>0) idx.00 <- idx.00[-unique(c(rm.1i,rm.2i))]
  
  if(length(idx.00)>0L){
    for (j in idx.00){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      prob<-ff[(idx.h+1):n1, (idx.v+1):n2]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))

      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl4.01<-rbind(pl4.01, nfyz12/dfyz12)
      }
    }
    
    pl4.sum<-colSums(pl4.01)/n
  } else{ pl4.sum <- 0 }
  
  pl1.sum+pl2.sum+pl3.sum+pl4.sum
  
}

pslogL.LM.p3<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  sigma.0<-parms$sigma.0
  err.y<-data$Y-theta%*%t(XZ)
  
  if(length(idx.111)>0L){
    pl.111.sum<-err.y[idx.111]%*%XZ[idx.111,]/n
  } else{ pl.111.sum <- 0 }
  
  pl.011 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.011,1]==n1)
  if( length(rm.i)>0) idx.011 <- idx.011[-rm.i]
  if(length(idx.011)>0L){
    for (j in idx.011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[(idx.h+1):n1, idx.v, idx.l]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j], data$z3st[j])
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.011<-rbind(pl.011, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z1.hat<-NULL
    }
    pl.011.sum<-colSums(pl.011)/n
    
  } else{ pl.011.sum <- 0 }
  
  pl.101 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.101, 2]==n2)
  if( length(rm.i)>0) idx.101 <- idx.101[-rm.i]
  if(length(idx.101)>0L){
    for (j in idx.101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[idx.h, (idx.v+1):n2, idx.l]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], z2.hat, data$z3st[j])
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.101<-rbind(pl.101, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z2.hat<-NULL
    }
    pl.101.sum<-colSums(pl.101)/n
    
  } else{ pl.101.sum <- 0 }
  
  pl.110 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.110, 3]==n3)
  if( length(rm.i)>0) idx.110 <- idx.110[-rm.i]
  if(length(idx.110)>0L){
    for (j in idx.110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[idx.h, idx.v, (idx.l+1):n3]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      ppq <- length(z3.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], z3.hat)
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.110<-rbind(pl.110, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z3.hat<-NULL
    }
    pl.110.sum<-colSums(pl.110)/n
    
  } else{ pl.110.sum <- 0 }
  
  pl.001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.001,1]==n1)
  rm.2i <- which(ord.res[idx.001,2]==n2)
  if( (length(rm.1i)+length(rm.2i))>0) idx.001 <- idx.001[-unique(c(rm.1i,rm.2i))]
  if(length(idx.001)>0L){
    for (j in idx.001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[(idx.h+1):n1, (idx.v+1):n2, idx.l]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1), data$z3st[j])
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.001<-rbind(pl.001, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
    }
    
    pl.001.sum<-colSums(pl.001)/n
  } else{ pl.001.sum <- 0 }
  
  pl.010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.010,1]==n1)
  rm.2i <- which(ord.res[idx.010,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.010 <- idx.010[-unique(c(rm.1i,rm.2i))]
  if(length(idx.010)>0L){
    for (j in idx.010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[(idx.h+1):n1, idx.v, (idx.l+1):n3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1))
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.010<-rbind(pl.010, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.010.sum<-colSums(pl.010)/n
  } else{ pl.010.sum <- 0 }
  
  pl.100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.100,2]==n2)
  rm.2i <- which(ord.res[idx.100,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.100 <- idx.100[-unique(c(rm.1i,rm.2i))]
  if(length(idx.100)>0L){
    for (j in idx.100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[idx.h, (idx.v+1):n2, (idx.l+1):n3]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1))
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.100<-rbind(pl.100, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.100.sum<-colSums(pl.100)/n
  } else{ pl.100.sum <- 0 }
  
  pl.000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.000,1]==n1)
  rm.2i <- which(ord.res[idx.000,2]==n2)
  rm.3i <- which(ord.res[idx.000,3]==n3)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.000 <- idx.000[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.000)>0L){
    for (j in idx.000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z3.hat)
      
      nt3 <- (idx.l+1):n3
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f123[(idx.h+1):n1, (idx.v+1):n2, k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k])
        
        # x.s <- outer(rep(1, p1*p2*p3), X[j,])
        # z123 <- cbind(matrix(rep(rbind(rep(z1.hat, p2), rep(z2.hat, each=p1)),p3) , p1*p2*p3, 2, byrow=T), rep(z3.hat, each=p1*p2))
        # XZ123 <- cbind(x.s, z123)
        
        Prob<-c(prob)
        err.yz123 <- c(data$Y[j]-theta%*%t(XZ123))
        fyz123 <- c(dnorm(err.yz123, sd=sigma.0))
        
        #nfyz123 <- c(fyz123*Prob*err.yz123)%*%XZ123
        #dfyz123 <- sum(fyz123*Prob)
        nfyz123 <- nfyz123 + c(fyz123*Prob*err.yz123)%*%XZ123
        dfyz123 <- dfyz123 + sum(fyz123*Prob)
        #nfyz123 <- rbind(nfyz123, c(fyz123*Prob*err.yz123)%*%XZ123)
        #dfyz123 <- c(dfyz123, sum(fyz123*Prob))
      }
      
      if(dfyz123!=0){
        pl.000<-rbind(pl.000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.000.sum<-colSums(pl.000)/n
    
  } else{ pl.000.sum <- 0 }
  
  pl.111.sum+pl.011.sum+pl.101.sum+pl.110.sum+pl.001.sum+pl.010.sum+pl.100.sum+pl.000.sum
  
}

pslogL.BN.p3<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(1/(1+exp(-ymu)))
  
  if(length(idx.111)>0L){
    pl.111.sum<-err.y[idx.111]%*%XZ[idx.111,]/n
  } else{ pl.111.sum <- 0 }
  
  pl.011 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.011,1]==n1)
  if( length(rm.i)>0) idx.011 <- idx.011[-rm.i]
  if(length(idx.011)>0L){
    for (j in idx.011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[(idx.h+1):n1, idx.v, idx.l]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j], data$z3st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.011<-rbind(pl.011, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z1.hat<-NULL
    }
    pl.011.sum<-colSums(pl.011)/n
    
  } else{ pl.011.sum <- 0 }
  
  pl.101 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.101, 2]==n2)
  if( length(rm.i)>0) idx.101 <- idx.101[-rm.i]
  if(length(idx.101)>0L){
    for (j in idx.101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[idx.h, (idx.v+1):n2, idx.l]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], z2.hat, data$z3st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.101<-rbind(pl.101, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z2.hat<-NULL
    }
    pl.101.sum<-colSums(pl.101)/n
    
  } else{ pl.101.sum <- 0 }
  
  pl.110 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.110, 3]==n3)
  if( length(rm.i)>0) idx.110 <- idx.110[-rm.i]
  if(length(idx.110)>0L){
    for (j in idx.110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[idx.h, idx.v, (idx.l+1):n3]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      ppq <- length(z3.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], z3.hat)
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.110<-rbind(pl.110, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z3.hat<-NULL
    }
    pl.110.sum<-colSums(pl.110)/n
    
  } else{ pl.110.sum <- 0 }
  
  pl.001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.001,1]==n1)
  rm.2i <- which(ord.res[idx.001,2]==n2)
  if( (length(rm.1i)+length(rm.2i))>0) idx.001 <- idx.001[-unique(c(rm.1i,rm.2i))]
  if(length(idx.001)>0L){
    for (j in idx.001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[(idx.h+1):n1, (idx.v+1):n2, idx.l]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1), data$z3st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.001<-rbind(pl.001, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
    }
    
    pl.001.sum<-colSums(pl.001)/n
  } else{ pl.001.sum <- 0 }
  
  pl.010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.010,1]==n1)
  rm.2i <- which(ord.res[idx.010,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.010 <- idx.010[-unique(c(rm.1i,rm.2i))]
  if(length(idx.010)>0L){
    for (j in idx.010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[(idx.h+1):n1, idx.v, (idx.l+1):n3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.010<-rbind(pl.010, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.010.sum<-colSums(pl.010)/n
  } else{ pl.010.sum <- 0 }
  
  pl.100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.100,2]==n2)
  rm.2i <- which(ord.res[idx.100,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.100 <- idx.100[-unique(c(rm.1i,rm.2i))]
  if(length(idx.100)>0L){
    for (j in idx.100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[idx.h, (idx.v+1):n2, (idx.l+1):n3]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.100<-rbind(pl.100, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.100.sum<-colSums(pl.100)/n
  } else{ pl.100.sum <- 0 }
  
  pl.000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.000,1]==n1)
  rm.2i <- which(ord.res[idx.000,2]==n2)
  rm.3i <- which(ord.res[idx.000,3]==n3)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.000 <- idx.000[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.000)>0L){
    for (j in idx.000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z3.hat)
      
      nt3 <- (idx.l+1):n3
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f123[(idx.h+1):n1, (idx.v+1):n2, k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k])
        
        # x.s <- outer(rep(1, p1*p2*p3), X[j,])
        # z123 <- cbind(matrix(rep(rbind(rep(z1.hat, p2), rep(z2.hat, each=p1)),p3) , p1*p2*p3, 2, byrow=T), rep(z3.hat, each=p1*p2))
        # XZ123 <- cbind(x.s, z123)
        
        Prob<-c(prob)
        y.mu <- theta%*%t(XZ123)
        y.pr <- c(1/(1+exp(-y.mu)))
        err.yz123<-c(data$Y[j]- y.pr)
        
        fyz123 <- c(dbinom(data$Y[j],1,prob=y.pr))

        #nfyz123 <- c(fyz123*Prob*err.yz123)%*%XZ123
        #dfyz123 <- sum(fyz123*Prob)
        nfyz123 <- nfyz123 + c(fyz123*Prob*err.yz123)%*%XZ123
        dfyz123 <- dfyz123 + sum(fyz123*Prob)
        #nfyz123 <- rbind(nfyz123, c(fyz123*Prob*err.yz123)%*%XZ123)
        #dfyz123 <- c(dfyz123, sum(fyz123*Prob))
      }
      
      if(dfyz123!=0){
        pl.000<-rbind(pl.000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.000.sum<-colSums(pl.000)/n
    
  } else{ pl.000.sum <- 0 }
  
  pl.111.sum+pl.011.sum+pl.101.sum+pl.110.sum+pl.001.sum+pl.010.sum+pl.100.sum+pl.000.sum
  
}

pslogL.PO.p3<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(exp(ymu))
  
  if(length(idx.111)>0L){
    pl.111.sum<-err.y[idx.111]%*%XZ[idx.111,]/n
  } else{ pl.111.sum <- 0 }
  
  pl.011 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.011,1]==n1)
  if( length(rm.i)>0) idx.011 <- idx.011[-rm.i]
  if(length(idx.011)>0L){
    for (j in idx.011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[(idx.h+1):n1, idx.v, idx.l]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j], data$z3st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.011<-rbind(pl.011, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z1.hat<-NULL
    }
    pl.011.sum<-colSums(pl.011)/n
    
  } else{ pl.011.sum <- 0 }
  
  pl.101 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.101, 2]==n2)
  if( length(rm.i)>0) idx.101 <- idx.101[-rm.i]
  if(length(idx.101)>0L){
    for (j in idx.101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[idx.h, (idx.v+1):n2, idx.l]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], z2.hat, data$z3st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.101<-rbind(pl.101, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z2.hat<-NULL
    }
    pl.101.sum<-colSums(pl.101)/n
    
  } else{ pl.101.sum <- 0 }
  
  pl.110 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.110, 3]==n3)
  if( length(rm.i)>0) idx.110 <- idx.110[-rm.i]
  if(length(idx.110)>0L){
    for (j in idx.110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      prob<-f123[idx.h, idx.v, (idx.l+1):n3]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      ppq <- length(z3.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], z3.hat)
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.110<-rbind(pl.110, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z3.hat<-NULL
    }
    pl.110.sum<-colSums(pl.110)/n
    
  } else{ pl.110.sum <- 0 }
  
  pl.001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.001,1]==n1)
  rm.2i <- which(ord.res[idx.001,2]==n2)
  if( (length(rm.1i)+length(rm.2i))>0) idx.001 <- idx.001[-unique(c(rm.1i,rm.2i))]
  if(length(idx.001)>0L){
    for (j in idx.001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[(idx.h+1):n1, (idx.v+1):n2, idx.l]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1), data$z3st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.001<-rbind(pl.001, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
    }
    
    pl.001.sum<-colSums(pl.001)/n
  } else{ pl.001.sum <- 0 }
  
  pl.010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.010,1]==n1)
  rm.2i <- which(ord.res[idx.010,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.010 <- idx.010[-unique(c(rm.1i,rm.2i))]
  if(length(idx.010)>0L){
    for (j in idx.010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[(idx.h+1):n1, idx.v, (idx.l+1):n3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.010<-rbind(pl.010, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.010.sum<-colSums(pl.010)/n
  } else{ pl.010.sum <- 0 }
  
  pl.100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.100,2]==n2)
  rm.2i <- which(ord.res[idx.100,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.100 <- idx.100[-unique(c(rm.1i,rm.2i))]
  if(length(idx.100)>0L){
    for (j in idx.100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      prob<-f123[idx.h, (idx.v+1):n2, (idx.l+1):n3]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.100<-rbind(pl.100, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.100.sum<-colSums(pl.100)/n
  } else{ pl.100.sum <- 0 }
  
  pl.000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.000,1]==n1)
  rm.2i <- which(ord.res[idx.000,2]==n2)
  rm.3i <- which(ord.res[idx.000,3]==n3)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.000 <- idx.000[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.000)>0L){
    for (j in idx.000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z3.hat)
      
      nt3 <- (idx.l+1):n3
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f123[(idx.h+1):n1, (idx.v+1):n2, k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k])
        
        # x.s <- outer(rep(1, p1*p2*p3), X[j,])
        # z123 <- cbind(matrix(rep(rbind(rep(z1.hat, p2), rep(z2.hat, each=p1)),p3) , p1*p2*p3, 2, byrow=T), rep(z3.hat, each=p1*p2))
        # XZ123 <- cbind(x.s, z123)
        
        Prob<-c(prob)
        y.mu <- theta%*%t(XZ123)
        y.pr <- c(exp(y.mu))
        err.yz123<-c(data$Y[j]- y.pr)
        
        fyz123 <- c(dpois(data$Y[j], lambda = y.pr))
        
        #nfyz123 <- c(fyz123*Prob*err.yz123)%*%XZ123
        #dfyz123 <- sum(fyz123*Prob)
        nfyz123 <- nfyz123 + c(fyz123*Prob*err.yz123)%*%XZ123
        dfyz123 <- dfyz123 + sum(fyz123*Prob)
        #nfyz123 <- rbind(nfyz123, c(fyz123*Prob*err.yz123)%*%XZ123)
        #dfyz123 <- c(dfyz123, sum(fyz123*Prob))
      }
      
      if(dfyz123!=0){
        pl.000<-rbind(pl.000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.000.sum<-colSums(pl.000)/n
    
  } else{ pl.000.sum <- 0 }
  
  pl.111.sum+pl.011.sum+pl.101.sum+pl.110.sum+pl.001.sum+pl.010.sum+pl.100.sum+pl.000.sum
  
}


# outer for two vectors given each subject
ind.outer <- function(x, y) {
  out<-array(0,c(ncol(x),dim(x)[1],dim(y)[1]))
  for(ss in 1:ncol(x)){
    out[ss,,]<-outer(x[,ss],y[,ss],"*")
  }
  return(out)
}

# outer for more than two vectors given each subject and eta(z(sj))>0
subj.outer <- function(ss, ind.use, prob.or, n.v, cenp, ind.k, subn) {
  
  plist <- list()
  nn<-which(ind.use[[paste0("Z",  ind.k[1])]][,subn[ss]]>0)
  m1<-ind.use[[paste0("Z",  ind.k[1])]][nn,subn[ss]]*prob.or[[paste0("p",  ind.k[1])]][nn]
  nn<-which(ind.use[[paste0("Z",  ind.k[2])]][,subn[ss]]>0)
  m2<-ind.use[[paste0("Z",  ind.k[2])]][nn,subn[ss]]*prob.or[[paste0("p",  ind.k[2])]][nn]
  plist[[1]] <- m1
  plist[[2]] <- m2
  len<-c(length(m1),length(m2))
  
  m1<-m2<-NULL
  for( us in 3:cenp){
    nn<-which(ind.use[[paste0("Z",  ind.k[us])]][,subn[ss]]>0)
    m1<-ind.use[[paste0("Z",  ind.k[us])]][nn,subn[ss]]*prob.or[[paste0("p",  ind.k[us])]][nn]
    plist[[us]] <- m1
    len<-c(len, length(m1))
  }
  
  pmat <- matrix(0, nrow=max(len), ncol=cenp)
  for (i in 1:cenp) {
    pmat[1:len[i],i] <- plist[[i]]
  }
  return(list(pmat=pmat, pdim=len))
}


# the dev. of log-pslikelihood which is the optimization function #
pslogL.LM.marg.mc<-function(theta, parms){
  # the parms includes:
  # data
  # XZ (1, X, Z)
  # miss.indx
  # est.z.r
  # res. or
  # prob.or
  # n.v
  # ord.res
  # ind.use
  # sigma.0
  # seedmc
  # all above will be obtained in the "main.func.LM.marg.mc" function
  
  data<-parms$data # matrix
  XZ<-as.matrix(parms$data$XZ) # matrix
  miss.indx<-parms$miss.indx # list ($miss1, $miss2, ... $p1, $p12, ...) of
  # vectors
  est.z.r<-parms$est.z.r # list ($Z1, $Z2, ...) of matrices
  res.or<-parms$res.or # list ($res1, $res2, ...) of vectors
  prob.or<-parms$prob.or # list ($p1, $p2, ...) of vectors
  n.v<-parms$n.v # vector
  ord.res<-parms$ord.res # matrix
  ind.use<-parms$ind.use # list ($Z1, $Z2, ...) of matrices
  sigma.0<-parms$sigma.0 # scalar
  
  err.y<-as.numeric(data$Y-theta%*%t(XZ))
  p<-dim(XZ)[2]
  pp<-dim(ord.res)[2]
  n<-dim(XZ)[1]
  qq<-p-pp
  X<-XZ[,c(1:qq)]
  
  ### complete-cases
  if(length(miss.indx$miss0)>0L){
    pl.cc<-err.y[miss.indx$miss0]%*%XZ[miss.indx$miss0,]/n
  } else{ pl.cc <- 0 }
  
  ### for one censoring
  pl.miss1 <- matrix(rep(0,p),1)
  cenp<-1
  
  cb<-combn(pp, cenp) # all combinations
  ncb<-dim(cb)[2]     # number of combinations
  for( jj in 1:ncb){
    ind.k<-cb[,jj]
    aa<-paste0(ind.k, collapse = "")
    subn<-miss.indx[[paste0("p",aa)]]  #subject in this combination
    
    if(length(subn)>0){
      rm.subj<-which(ord.res[subn, ind.k]==n.v[ind.k])
      if(length(rm.subj)>0)  subn<-subn[-rm.subj]
      if(length(subn)>0){
        
        # prob=eta(ds)
        prob<-matrix(ind.use[[paste0("Z",  ind.k)]][,subn]*prob.or[[paste0("p",  ind.k)]],,length(subn))
        # est.z(s)
        zk.hat <- matrix(est.z.r[[paste0("Z", ind.k)]][,subn],,length(subn))
        ppq <- dim(zk.hat)[1] # number of sum points
        
        # error term (y-hat.y)
        err.yz11 <- c(data$Y[subn] - theta[-c(qq+ind.k)]%*%t(matrix(XZ[subn, -c(qq+ind.k)],length(subn))) )
        err.yz1 <- matrix(err.yz11,ppq,length(subn),byrow=T)-theta[c(qq+ind.k)]*zk.hat
        
        fyz1 <- dnorm(err.yz1, sd=sigma.0) # f(y|x,z,z(s))
        
        if(sum(fyz1)>0){
          nfyz10 <- fyz1*prob*err.yz1        # f(y|x,z,z(s))*eta(s)*err*1
          nfyz11 <- fyz1*prob*err.yz1*zk.hat # f(y|x,z,z(s))*eta(s)*err*z(s)
          
          # nfyz1 for intercept and z(s) = numerator = sum(f(y|x,z,z(s))*eta(s)*err*(1, z(s))^T)
          nfyz1 <- cbind(colSums(nfyz10),colSums(nfyz11))
          # dfyz1 for intercept and z(s) = denominator = sum(f(y|x,z,z(s))*eta(s))
          dfyz1 <- colSums(fyz1*prob)
          
          ps<-which(dfyz1>0)
          if(length(ps)>0){
            ap <- matrix((nfyz1/dfyz1)[ps,],length(subn[ps]))
            
            # make nfyz1/dfyz1 for all design matrix (1, x1...xq, zx,....,z(s),..): dxz
            dxz <- matrix(XZ[subn[ps], ],length(subn[ps]))
            dxz[, -c(qq+ind.k)]<-matrix(XZ[subn[ps], -c(qq+ind.k)],length(subn[ps]))*ap[,1]
            dxz[, c(qq+ind.k)]<-ap[,2]
            
            pl.miss1<-rbind(pl.miss1, dxz) #dev.pslog-likelihood for subjects with one censoring covariat z(s)
            
          }
          
          prob<-zk.hat<-err.y1<-fyz1<-nfyz10<-nfyz11<-NULL
          
        }
      }
      
    }
    
  }
  
  if(dim(pl.miss1)[1]>1){
    pl.m1 <- colSums(pl.miss1)/n
  }else{ pl.m1<-0}
  
  ### for subejcts with two censorings covariates
  pl.miss2 <- matrix(rep(0,p),1)
  cenp<-2
  
  cb<-combn(pp, cenp)
  ncb<-dim(cb)[2]
  for( jj in 1:ncb){
    ind.k<-cb[,jj]
    aa<-paste0(ind.k, collapse = "")
    subn<-miss.indx[[paste0("p",aa)]]
    
    if(length(subn)>0){
      rm.subj<-c()
      for( us in 1:cenp){
        rm.subj<-c(rm.subj,which(ord.res[subn, ind.k[us]]==n.v[ind.k[us]]))
      }
      if(length(rm.subj)>0)  subn<-subn[-unique(rm.subj)]
      
      if(length(subn)>0){
        prob<-ind.outer(x=matrix(ind.use[[paste0("Z",  ind.k[1])]][,subn]*prob.or[[paste0("p",  ind.k[1])]],,length(subn)),
                        y=matrix(ind.use[[paste0("Z",  ind.k[2])]][,subn]*prob.or[[paste0("p",  ind.k[2])]],,length(subn)))
        
        zk.hat1 <- aperm(array(matrix(est.z.r[[paste0("Z", ind.k[1])]][,subn],,length(subn)), c(n.v[ind.k[1]],length(subn),n.v[ind.k[2]] )), c(2,1,3))
        zk.hat2 <- aperm(array(matrix(est.z.r[[paste0("Z", ind.k[2])]][,subn],,length(subn)), c(n.v[ind.k[2]],length(subn),n.v[ind.k[1]] )), c(2,3,1))
        
        err.yz11 <- c(data$Y[subn] - theta[-c(qq+ind.k)]%*%t(matrix(XZ[subn, -c(qq+ind.k)],length(subn))  ))
        
        err.yz1 <- err.yz11-theta[c(qq+ind.k[1])]*zk.hat1-theta[c(qq+ind.k[2])]*zk.hat2
        
        fyz1 <- dnorm(err.yz1, sd=sigma.0)
        
        if(sum(fyz1)>0){
          nfyz10 <- fyz1*prob*err.yz1         # f(y|x,z,z(s))*eta(s)*err*1
          nfyz11 <- fyz1*prob*err.yz1*zk.hat1 # f(y|x,z,z(s))*eta(s)*err*z(s1)
          nfyz12 <- fyz1*prob*err.yz1*zk.hat2 # f(y|x,z,z(s))*eta(s)*err*z(s2)
          
          nfyz1 <- cbind(apply(nfyz10,1,sum),apply(nfyz11,1,sum),apply(nfyz12,1,sum))
          dfyz1 <- apply(fyz1*prob,1,sum)
          
          ps<-which(dfyz1>0)
          
          if(length(ps)>0){
            ap <- matrix((nfyz1/dfyz1)[ps,],length(subn[ps]))
            
            # make nfyz1/dfyz1 for all design matrix (1, x1...xq, zx,..z(s1)..,z(s2),..) :dxz
            dxz <- matrix(XZ[subn[ps], ],length(subn[ps]))
            dxz[, c(qq+ind.k)]<-ap[,-1]
            dxz[, -c(qq+ind.k)]<-matrix(XZ[subn[ps], -c(qq+ind.k)],length(subn[ps]))*ap[,1]
            
            pl.miss2<-rbind(pl.miss2, dxz) #dev.pslog-likelihood for subjects with two censoring covariats z(s1), z(s2)
            
          }
          
          prob<-zk.hat1<-zk.hat2<-err.y1<-fyz1<-nfyz10<-nfyz11<-nfyz12<-NULL
          
        }
        
      }
      
    }
    
    
  }
  
  if(dim(pl.miss2)[1]>1){
    pl.m2 <- colSums(pl.miss2)/n
  }else{ pl.m2<-0}
  
  #### subjects with more than 2 covariates subject to LOD
  pl.missmore2 <- matrix(rep(0,p),1)
  
  for(cenp in 3:pp){ # number of censoring covariates (cenp) from 3 up to pp
    cb<-combn(pp, cenp) # all combinations under cenp
    ncb<-dim(cb)[2]     # number of combinations
    
    for( jj in 1:ncb){
      ind.k<-cb[,jj]
      aa<-paste0(ind.k, collapse = "")
      subn<-miss.indx[[paste0("p",aa)]]
      
      if(length(subn)>0){
        rm.subj<-c()
        for( us in 1:cenp){
          rm.subj<-c(rm.subj,which(ord.res[subn, ind.k[us]]==n.v[ind.k[us]]))
        }
        if(length(rm.subj)>0)  subn<-subn[-unique(rm.subj)]
        
        if(length(subn)>0){
          err.yz11.a <- c(data$Y[subn] - theta[-c(qq+ind.k)]%*%t(matrix(XZ[subn, -c(qq+ind.k)],length(subn))  ))
          
          # due to the memory issue,
          # we here calcuate the dev.log-likelihood for each subject, instead of a vector of all subjects
          for( ss in 1:length(subn)){
            pout<-subj.outer(ss, ind.use, prob.or, n.v, cenp, ind.k, subn)
            zk <- calc_zk(cenp, ss, subn, ind.k, ord.res, est.z.r,
                          pout$pdim)
            if (cumprod(pout$pdim)[length(pout$pdim)]>1000000) {
              
              # given a fixed seed to make sure the drawn samples are the same.
              set.seed(parms$seedmc)
              
              fyztmp <- calc_nfyz1_mc(err.yz11.a[ss], qq, ind.k, theta, zk,
                                      pout$pmat, sigma.0,
                                      pout$pdim, 100000)
            }
            else {
              fyztmp <- calc_nfyz1(err.yz11.a[ss], qq, ind.k, theta, zk,
                                   pout$pmat, pout$pdim, sigma.0)
            }
            
            if( sum(fyztmp)>0){
              nfyz1 <- fyztmp[1:(length(fyztmp)-1)]
              dfyz1 <- fyztmp[length(fyztmp)]
              
              if(dfyz1>0){
                ap<-nfyz1/dfyz1
                dxz<-XZ[subn[ss], ]
                dxz[-c(qq+ind.k)]<-XZ[subn[ss], -c(qq+ind.k)]*ap[1]
                dxz[c(qq+ind.k)]<-ap[-1]
                pl.missmore2<-rbind(pl.missmore2, dxz) # dev.log-likelihood for all subjects with more than 2 censoring covariates
              }
            }
            
          }
        }
        
      }
    }
  }
  
  if(dim(pl.missmore2)[1]>1){
    pl.m3 <- colSums(pl.missmore2)/n
  }else{ pl.m3<-0}
  
  pl.cc+pl.m1+pl.m2+pl.m3
}
