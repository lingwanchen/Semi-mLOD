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

  if(length(data$idx.11)>0L){
    pl1.sum<-err.y[data$idx.11]%*%XZ[data$idx.11,]/n
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

pslogL.LM.marg.p4<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  sigma.0<-parms$sigma.0
  err.y<-data$Y-theta%*%t(XZ)
  
  if(length(idx.1111)>0L){
    pl.1111.sum<-err.y[idx.1111]%*%XZ[idx.1111,]/n
  } else{ pl.1111.sum <- 0 }
  
  pl.0111 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.0111,1]==n1)
  if( length(rm.i)>0) idx.0111 <- idx.0111[-rm.i]
  if(length(idx.0111)>0L){
    for (j in idx.0111){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f1[(idx.h+1):n1]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j], data$z3st[j], data$z4st[j])
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.0111<-rbind(pl.0111, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z1.hat<-NULL
    }
    pl.0111.sum<-colSums(pl.0111)/n
    
  } else{ pl.0111.sum <- 0 }
  
  pl.1011 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1011, 2]==n2)
  if( length(rm.i)>0) idx.1011 <- idx.1011[-rm.i]
  if(length(idx.1011)>0L){
    for (j in idx.1011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f2[(idx.v+1):n2]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], z2.hat, data$z3st[j], data$z4st[j])
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1011<-rbind(pl.1011, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z2.hat<-NULL
    }
    pl.1011.sum<-colSums(pl.1011)/n
    
  } else{ pl.1011.sum <- 0 }
  
  pl.1101 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1101, 3]==n3)
  if( length(rm.i)>0) idx.1101 <- idx.1101[-rm.i]
  if(length(idx.1101)>0L){
    for (j in idx.1101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f3[(idx.l+1):n3]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      ppq <- length(z3.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], z3.hat, data$z4st[j])
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1101<-rbind(pl.1101, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z3.hat<-NULL
    }
    pl.1101.sum<-colSums(pl.1101)/n
    
  } else{ pl.1101.sum <- 0 }
  
  pl.1110 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1110, 4]==n4)
  if( length(rm.i)>0) idx.1110 <- idx.1110[-rm.i]
  if(length(idx.1110)>0L){
    for (j in idx.1110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f4[(idx.k+1):n4]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      ppq <- length(z4.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], data$z3st[j], z4.hat)
      
      err.yz1 <- c(data$Y[j] - theta%*%t(XZ1))
      fyz1 <- c(dnorm(err.yz1, sd=sigma.0))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1110<-rbind(pl.1110, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z4.hat<-NULL
    }
    pl.1110.sum<-colSums(pl.1110)/n
    
  } else{ pl.1110.sum <- 0 }
  
  
  pl.0011 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0011,1]==n1)
  rm.2i <- which(ord.res[idx.0011,2]==n2)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0011 <- idx.0011[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0011)>0L){
    for (j in idx.0011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f12[(idx.h+1):n1, (idx.v+1):n2]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1), data$z3st[j], data$z4st[j])
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0011<-rbind(pl.0011, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
    }
    
    pl.0011.sum<-colSums(pl.0011)/n
  } else{ pl.0011.sum <- 0 }
  
  pl.0101 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0101,1]==n1)
  rm.2i <- which(ord.res[idx.0101,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0101 <- idx.0101[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0101)>0L){
    for (j in idx.0101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f13[(idx.h+1):n1,(idx.l+1):n3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1), data$z4st[j])
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0101<-rbind(pl.0101, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.0101.sum<-colSums(pl.0101)/n
  } else{ pl.0101.sum <- 0 }
  
  pl.0110 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0110,1]==n1)
  rm.2i <- which(ord.res[idx.0110,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0110 <- idx.0110[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0110)>0L){
    for (j in idx.0110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f14[(idx.h+1):n1, (idx.k+1):n4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], data$z3st[j], rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0110<-rbind(pl.0110, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0110.sum<-colSums(pl.0110)/n
  } else{ pl.0110.sum <- 0 }
  
  pl.1001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1001,2]==n2)
  rm.2i <- which(ord.res[idx.1001,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1001 <- idx.1001[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1001)>0L){
    for (j in idx.1001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f23[(idx.v+1):n2, (idx.l+1):n3]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1), data$z4st[j])
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1001<-rbind(pl.1001, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.1001.sum<-colSums(pl.1001)/n
  } else{ pl.1001.sum <- 0 }
  
  pl.1010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1010,2]==n2)
  rm.2i <- which(ord.res[idx.1010,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1010 <- idx.1010[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1010)>0L){
    for (j in idx.1010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f24[(idx.v+1):n2, (idx.k+1):n4]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), data$z3st[j], rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1010<-rbind(pl.1010, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1010.sum<-colSums(pl.1010)/n
  } else{ pl.1010.sum <- 0 }
  
  pl.1100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1100,3]==n3)
  rm.2i <- which(ord.res[idx.1100,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1100 <- idx.1100[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1100)>0L){
    for (j in idx.1100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f34[(idx.l+1):n3, (idx.k+1):n4]
      
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z3.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], data$z2st[j], rep(z3.hat, p2), rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      err.yz12 <- c(data$Y[j]-theta%*%t(XZ12))
      fyz12 <- c(dnorm(err.yz12, sd=sigma.0))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1100<-rbind(pl.1100, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1100.sum<-colSums(pl.1100)/n
  } else{ pl.1100.sum <- 0 }
  
  pl.0001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0001,1]==n1)
  rm.2i <- which(ord.res[idx.0001,2]==n2)
  rm.3i <- which(ord.res[idx.0001,3]==n3)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0001 <- idx.0001[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0001)>0L){
    for (j in idx.0001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
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
        prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f3[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k], data$z4st[j])
        
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
        pl.0001<-rbind(pl.0001, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.0001.sum<-colSums(pl.0001)/n
    
  } else{ pl.0001.sum <- 0 }
  
  pl.0010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0010,1]==n1)
  rm.2i <- which(ord.res[idx.0010,2]==n2)
  rm.3i <- which(ord.res[idx.0010,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0010 <- idx.0010[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0010)>0L){
    for (j in idx.0010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), data$z3st[j], z4.hat[id.k])
        
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
        pl.0010<-rbind(pl.0010, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0010.sum<-colSums(pl.0010)/n
    
  } else{ pl.0010.sum <- 0 }
  
  pl.0100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0100,1]==n1)
  rm.2i <- which(ord.res[idx.0100,3]==n3)
  rm.3i <- which(ord.res[idx.0100,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0100 <- idx.0100[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0100)>0L){
    for (j in idx.0100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f13[(idx.h+1):n1, (idx.l+1):n3]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1), z4.hat[id.k])
        
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
        pl.0100<-rbind(pl.0100, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0100.sum<-colSums(pl.0100)/n
    
  } else{ pl.0100.sum <- 0 }
  
  pl.1000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1000,2]==n2)
  rm.2i <- which(ord.res[idx.1000,3]==n3)
  rm.3i <- which(ord.res[idx.1000,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.1000 <- idx.1000[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.1000)>0L){
    for (j in idx.1000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f23[(idx.v+1):n2, (idx.l+1):n3]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1), z4.hat[id.k])
        
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
        pl.1000<-rbind(pl.1000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1000.sum<-colSums(pl.1000)/n
    
  } else{ pl.1000.sum <- 0 }
  
  
  pl.0000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0000,1]==n1)
  rm.2i <- which(ord.res[idx.0000,2]==n2)
  rm.3i <- which(ord.res[idx.0000,3]==n3)
  rm.4i <- which(ord.res[idx.0000,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i)+length(rm.4i))>0) idx.0000 <- idx.0000[-unique(c(rm.1i,rm.2i,rm.3i,rm.4i))]
  if(length(idx.0000)>0L){
    for (j in idx.0000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z3.hat)
      p4 <- length(z4.hat)
      
      nt3 <- (idx.l+1):n3
      nt4 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.kk in 1:p4){
        
        kk <- nt4[id.kk]
        
        for( id.k in 1:p3){
          
          k <- nt3[id.k]
          prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f3[k]*f4[kk]
          
          x.s <- outer(rep(1,p1*p2), X[j,])
          XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k], z4.hat[id.kk])
          
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
      }
      
      if(dfyz123!=0){
        pl.0000<-rbind(pl.0000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0000.sum<-colSums(pl.0000)/n
    
  } else{ pl.0000.sum <- 0 }
  
  pl.1111.sum+pl.0111.sum+pl.1011.sum+pl.1101.sum+pl.1110.sum+
    pl.0011.sum+pl.0101.sum+pl.0110.sum+pl.1001.sum+pl.1010.sum+pl.1100.sum+
    pl.0001.sum+pl.0010.sum+pl.0100.sum+pl.1000.sum+pl.0000.sum
  
}

pslogL.BN.marg.p4<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(1/(1+exp(-ymu)))
  
  if(length(idx.1111)>0L){
    pl.1111.sum<-err.y[idx.1111]%*%XZ[idx.1111,]/n
  } else{ pl.1111.sum <- 0 }
  
  pl.0111 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.0111,1]==n1)
  if( length(rm.i)>0) idx.0111 <- idx.0111[-rm.i]
  if(length(idx.0111)>0L){
    for (j in idx.0111){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f1[(idx.h+1):n1]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j], data$z3st[j], data$z4st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.0111<-rbind(pl.0111, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z1.hat<-NULL
    }
    pl.0111.sum<-colSums(pl.0111)/n
    
  } else{ pl.0111.sum <- 0 }
  
  pl.1011 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1011, 2]==n2)
  if( length(rm.i)>0) idx.1011 <- idx.1011[-rm.i]
  if(length(idx.1011)>0L){
    for (j in idx.1011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f2[(idx.v+1):n2]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], z2.hat, data$z3st[j], data$z4st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1011<-rbind(pl.1011, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z2.hat<-NULL
    }
    pl.1011.sum<-colSums(pl.1011)/n
    
  } else{ pl.1011.sum <- 0 }
  
  pl.1101 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1101, 3]==n3)
  if( length(rm.i)>0) idx.1101 <- idx.1101[-rm.i]
  if(length(idx.1101)>0L){
    for (j in idx.1101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f3[(idx.l+1):n3]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      ppq <- length(z3.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], z3.hat, data$z4st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1101<-rbind(pl.1101, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z3.hat<-NULL
    }
    pl.1101.sum<-colSums(pl.1101)/n
    
  } else{ pl.1101.sum <- 0 }
  
  pl.1110 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1110, 4]==n4)
  if( length(rm.i)>0) idx.1110 <- idx.1110[-rm.i]
  if(length(idx.1110)>0L){
    for (j in idx.1110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f4[(idx.k+1):n4]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      ppq <- length(z4.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], data$z3st[j], z4.hat)
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1110<-rbind(pl.1110, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z4.hat<-NULL
    }
    pl.1110.sum<-colSums(pl.1110)/n
    
  } else{ pl.1110.sum <- 0 }
  
  
  pl.0011 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0011,1]==n1)
  rm.2i <- which(ord.res[idx.0011,2]==n2)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0011 <- idx.0011[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0011)>0L){
    for (j in idx.0011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f12[(idx.h+1):n1, (idx.v+1):n2]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1), data$z3st[j], data$z4st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0011<-rbind(pl.0011, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
    }
    
    pl.0011.sum<-colSums(pl.0011)/n
  } else{ pl.0011.sum <- 0 }
  
  pl.0101 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0101,1]==n1)
  rm.2i <- which(ord.res[idx.0101,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0101 <- idx.0101[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0101)>0L){
    for (j in idx.0101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f13[(idx.h+1):n1,(idx.l+1):n3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1), data$z4st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0101<-rbind(pl.0101, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.0101.sum<-colSums(pl.0101)/n
  } else{ pl.0101.sum <- 0 }
  
  pl.0110 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0110,1]==n1)
  rm.2i <- which(ord.res[idx.0110,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0110 <- idx.0110[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0110)>0L){
    for (j in idx.0110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f14[(idx.h+1):n1, (idx.k+1):n4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], data$z3st[j], rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0110<-rbind(pl.0110, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0110.sum<-colSums(pl.0110)/n
  } else{ pl.0110.sum <- 0 }
  
  pl.1001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1001,2]==n2)
  rm.2i <- which(ord.res[idx.1001,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1001 <- idx.1001[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1001)>0L){
    for (j in idx.1001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f23[(idx.v+1):n2, (idx.l+1):n3]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1), data$z4st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1001<-rbind(pl.1001, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.1001.sum<-colSums(pl.1001)/n
  } else{ pl.1001.sum <- 0 }
  
  pl.1010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1010,2]==n2)
  rm.2i <- which(ord.res[idx.1010,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1010 <- idx.1010[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1010)>0L){
    for (j in idx.1010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f24[(idx.v+1):n2, (idx.k+1):n4]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), data$z3st[j], rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1010<-rbind(pl.1010, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1010.sum<-colSums(pl.1010)/n
  } else{ pl.1010.sum <- 0 }
  
  pl.1100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1100,3]==n3)
  rm.2i <- which(ord.res[idx.1100,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1100 <- idx.1100[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1100)>0L){
    for (j in idx.1100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f34[(idx.l+1):n3, (idx.k+1):n4]
      
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z3.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], data$z2st[j], rep(z3.hat, p2), rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1100<-rbind(pl.1100, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1100.sum<-colSums(pl.1100)/n
  } else{ pl.1100.sum <- 0 }
  
  pl.0001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0001,1]==n1)
  rm.2i <- which(ord.res[idx.0001,2]==n2)
  rm.3i <- which(ord.res[idx.0001,3]==n3)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0001 <- idx.0001[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0001)>0L){
    for (j in idx.0001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
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
        prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f3[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k], data$z4st[j])
        
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
        pl.0001<-rbind(pl.0001, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.0001.sum<-colSums(pl.0001)/n
    
  } else{ pl.0001.sum <- 0 }
  
  pl.0010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0010,1]==n1)
  rm.2i <- which(ord.res[idx.0010,2]==n2)
  rm.3i <- which(ord.res[idx.0010,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0010 <- idx.0010[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0010)>0L){
    for (j in idx.0010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), data$z3st[j], z4.hat[id.k])
        
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
        pl.0010<-rbind(pl.0010, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0010.sum<-colSums(pl.0010)/n
    
  } else{ pl.0010.sum <- 0 }
  
  pl.0100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0100,1]==n1)
  rm.2i <- which(ord.res[idx.0100,3]==n3)
  rm.3i <- which(ord.res[idx.0100,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0100 <- idx.0100[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0100)>0L){
    for (j in idx.0100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f13[(idx.h+1):n1, (idx.l+1):n3]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1), z4.hat[id.k])
        
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
        pl.0100<-rbind(pl.0100, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0100.sum<-colSums(pl.0100)/n
    
  } else{ pl.0100.sum <- 0 }
  
  pl.1000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1000,2]==n2)
  rm.2i <- which(ord.res[idx.1000,3]==n3)
  rm.3i <- which(ord.res[idx.1000,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.1000 <- idx.1000[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.1000)>0L){
    for (j in idx.1000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f23[(idx.v+1):n2, (idx.l+1):n3]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1), z4.hat[id.k])
        
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
        pl.1000<-rbind(pl.1000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1000.sum<-colSums(pl.1000)/n
    
  } else{ pl.1000.sum <- 0 }
  
  
  pl.0000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0000,1]==n1)
  rm.2i <- which(ord.res[idx.0000,2]==n2)
  rm.3i <- which(ord.res[idx.0000,3]==n3)
  rm.4i <- which(ord.res[idx.0000,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i)+length(rm.4i))>0) idx.0000 <- idx.0000[-unique(c(rm.1i,rm.2i,rm.3i,rm.4i))]
  if(length(idx.0000)>0L){
    for (j in idx.0000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z3.hat)
      p4 <- length(z4.hat)
      
      nt3 <- (idx.l+1):n3
      nt4 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.kk in 1:p4){
        
        kk <- nt4[id.kk]
        
        for( id.k in 1:p3){
          
          k <- nt3[id.k]
          prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f3[k]*f4[kk]
          
          x.s <- outer(rep(1,p1*p2), X[j,])
          XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k], z4.hat[id.kk])
          
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
      }
      
      if(dfyz123!=0){
        pl.0000<-rbind(pl.0000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0000.sum<-colSums(pl.0000)/n
    
  } else{ pl.0000.sum <- 0 }
  
  pl.1111.sum+pl.0111.sum+pl.1011.sum+pl.1101.sum+pl.1110.sum+
    pl.0011.sum+pl.0101.sum+pl.0110.sum+pl.1001.sum+pl.1010.sum+pl.1100.sum+
    pl.0001.sum+pl.0010.sum+pl.0100.sum+pl.1000.sum+pl.0000.sum
  
}

pslogL.PO.marg.p4<-function(theta, parms){
  
  data<-parms$data
  XZ<-data$XZ
  n<-parms$n
  p<-parms$p
  X<-data$XX
  ymu <- theta%*%t(XZ)
  err.y<-data$Y- c(exp(ymu))
  
  if(length(idx.1111)>0L){
    pl.1111.sum<-err.y[idx.1111]%*%XZ[idx.1111,]/n
  } else{ pl.1111.sum <- 0 }
  
  pl.0111 <- matrix(rep(0,p),1)
  rm.i <- which(ord.res[idx.0111,1]==n1)
  if( length(rm.i)>0) idx.0111 <- idx.0111[-rm.i]
  if(length(idx.0111)>0L){
    for (j in idx.0111){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f1[(idx.h+1):n1]
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      
      ppq <- length(z1.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, z1.hat, data$z2st[j], data$z3st[j], data$z4st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.0111<-rbind(pl.0111, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z1.hat<-NULL
    }
    pl.0111.sum<-colSums(pl.0111)/n
    
  } else{ pl.0111.sum <- 0 }
  
  pl.1011 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1011, 2]==n2)
  if( length(rm.i)>0) idx.1011 <- idx.1011[-rm.i]
  if(length(idx.1011)>0L){
    for (j in idx.1011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f2[(idx.v+1):n2]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      ppq <- length(z2.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], z2.hat, data$z3st[j], data$z4st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1011<-rbind(pl.1011, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z2.hat<-NULL
    }
    pl.1011.sum<-colSums(pl.1011)/n
    
  } else{ pl.1011.sum <- 0 }
  
  pl.1101 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1101, 3]==n3)
  if( length(rm.i)>0) idx.1101 <- idx.1101[-rm.i]
  if(length(idx.1101)>0L){
    for (j in idx.1101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f3[(idx.l+1):n3]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      ppq <- length(z3.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], z3.hat, data$z4st[j])
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1101<-rbind(pl.1101, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z3.hat<-NULL
    }
    pl.1101.sum<-colSums(pl.1101)/n
    
  } else{ pl.1101.sum <- 0 }
  
  pl.1110 <- matrix(rep(0,p), 1)
  rm.i <- which(ord.res[idx.1110, 4]==n4)
  if( length(rm.i)>0) idx.1110 <- idx.1110[-rm.i]
  if(length(idx.1110)>0L){
    for (j in idx.1110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      prob<-f4[(idx.k+1):n4]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      ppq <- length(z4.hat)
      x.s <- outer(rep(1,ppq), X[j,]) # X includes int.
      XZ1 <- cbind(x.s, data$z1st[j], data$z2st[j], data$z3st[j], z4.hat)
      
      y.mu <- theta%*%t(XZ1)
      y.pr <- c(exp(y.mu))
      err.yz1<-c(data$Y[j]- y.pr)
      
      fyz1 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz1 <- c(fyz1*prob*err.yz1)%*%XZ1
      dfyz1 <- sum(fyz1*prob)
      if(dfyz1>0L){
        pl.1110<-rbind(pl.1110, nfyz1/dfyz1)
      }
      x.s<-NULL
      XZ1<-NULL
      prob<-NULL
      z4.hat<-NULL
    }
    pl.1110.sum<-colSums(pl.1110)/n
    
  } else{ pl.1110.sum <- 0 }
  
  
  pl.0011 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0011,1]==n1)
  rm.2i <- which(ord.res[idx.0011,2]==n2)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0011 <- idx.0011[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0011)>0L){
    for (j in idx.0011){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f12[(idx.h+1):n1, (idx.v+1):n2]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat,p2), rep(z2.hat,each=p1), data$z3st[j], data$z4st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0011<-rbind(pl.0011, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
    }
    
    pl.0011.sum<-colSums(pl.0011)/n
  } else{ pl.0011.sum <- 0 }
  
  pl.0101 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0101,1]==n1)
  rm.2i <- which(ord.res[idx.0101,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0101 <- idx.0101[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0101)>0L){
    for (j in idx.0101){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f13[(idx.h+1):n1,(idx.l+1):n3]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1), data$z4st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0101<-rbind(pl.0101, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.0101.sum<-colSums(pl.0101)/n
  } else{ pl.0101.sum <- 0 }
  
  pl.0110 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0110,1]==n1)
  rm.2i <- which(ord.res[idx.0110,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.0110 <- idx.0110[-unique(c(rm.1i,rm.2i))]
  if(length(idx.0110)>0L){
    for (j in idx.0110){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f14[(idx.h+1):n1, (idx.k+1):n4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], data$z3st[j], rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.0110<-rbind(pl.0110, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z1.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0110.sum<-colSums(pl.0110)/n
  } else{ pl.0110.sum <- 0 }
  
  pl.1001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1001,2]==n2)
  rm.2i <- which(ord.res[idx.1001,3]==n3)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1001 <- idx.1001[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1001)>0L){
    for (j in idx.1001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f23[(idx.v+1):n2, (idx.l+1):n3]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1), data$z4st[j])
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1001<-rbind(pl.1001, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.1001.sum<-colSums(pl.1001)/n
  } else{ pl.1001.sum <- 0 }
  
  pl.1010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1010,2]==n2)
  rm.2i <- which(ord.res[idx.1010,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1010 <- idx.1010[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1010)>0L){
    for (j in idx.1010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f24[(idx.v+1):n2, (idx.k+1):n4]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), data$z3st[j], rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(1/(1+exp(-y.mu)))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dbinom(data$Y[j],1,prob=y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1010<-rbind(pl.1010, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z2.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1010.sum<-colSums(pl.1010)/n
  } else{ pl.1010.sum <- 0 }
  
  pl.1100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1100,3]==n3)
  rm.2i <- which(ord.res[idx.1100,4]==n4)
  if( (length(rm.1i)+length(rm.2i))>0) idx.1100 <- idx.1100[-unique(c(rm.1i,rm.2i))]
  if(length(idx.1100)>0L){
    for (j in idx.1100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      prob<-f34[(idx.l+1):n3, (idx.k+1):n4]
      
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z3.hat)
      p2 <- length(z4.hat)
      
      x.s <- outer(rep(1,p1*p2), X[j,])
      XZ12 <- cbind(x.s, data$z1st[j], data$z2st[j], rep(z3.hat, p2), rep(z4.hat, each=p1))
      Prob<-c(prob)
      
      y.mu <- theta%*%t(XZ12)
      y.pr <- c(exp(y.mu))
      err.yz12<-c(data$Y[j]- y.pr)
      
      fyz12 <- c(dpois(data$Y[j], lambda = y.pr))
      nfyz12 <- c(fyz12*Prob*err.yz12)%*%XZ12
      #te <-fyz12*Prob
      dfyz12 <- sum(fyz12*Prob)
      if(dfyz12!=0){
        pl.1100<-rbind(pl.1100, nfyz12/dfyz12)
      }
      x.s<-NULL
      XZ12<-NULL
      prob<-NULL
      Prob<-err.yz12<-fyz12<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1100.sum<-colSums(pl.1100)/n
  } else{ pl.1100.sum <- 0 }
  
  pl.0001 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0001,1]==n1)
  rm.2i <- which(ord.res[idx.0001,2]==n2)
  rm.3i <- which(ord.res[idx.0001,3]==n3)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0001 <- idx.0001[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0001)>0L){
    for (j in idx.0001){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
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
        prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f3[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k], data$z4st[j])
        
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
        pl.0001<-rbind(pl.0001, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
    }
    
    pl.0001.sum<-colSums(pl.0001)/n
    
  } else{ pl.0001.sum <- 0 }
  
  pl.0010 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0010,1]==n1)
  rm.2i <- which(ord.res[idx.0010,2]==n2)
  rm.3i <- which(ord.res[idx.0010,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0010 <- idx.0010[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0010)>0L){
    for (j in idx.0010){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), data$z3st[j], z4.hat[id.k])
        
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
        pl.0010<-rbind(pl.0010, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0010.sum<-colSums(pl.0010)/n
    
  } else{ pl.0010.sum <- 0 }
  
  pl.0100 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0100,1]==n1)
  rm.2i <- which(ord.res[idx.0100,3]==n3)
  rm.3i <- which(ord.res[idx.0100,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.0100 <- idx.0100[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.0100)>0L){
    for (j in idx.0100){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z3.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f13[(idx.h+1):n1, (idx.l+1):n3]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, rep(z1.hat, p2), data$z2st[j], rep(z3.hat, each=p1), z4.hat[id.k])
        
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
        pl.0100<-rbind(pl.0100, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0100.sum<-colSums(pl.0100)/n
    
  } else{ pl.0100.sum <- 0 }
  
  pl.1000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.1000,2]==n2)
  rm.2i <- which(ord.res[idx.1000,3]==n3)
  rm.3i <- which(ord.res[idx.1000,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i))>0) idx.1000 <- idx.1000[-unique(c(rm.1i,rm.2i,rm.3i))]
  if(length(idx.1000)>0L){
    for (j in idx.1000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z2.hat)
      p2 <- length(z3.hat)
      p3 <- length(z4.hat)
      
      nt3 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.k in 1:p3){
        
        k <- nt3[id.k]
        prob<-f23[(idx.v+1):n2, (idx.l+1):n3]*f4[k]
        
        x.s <- outer(rep(1,p1*p2), X[j,])
        XZ123 <- cbind(x.s, data$z1st[j], rep(z2.hat, p2), rep(z3.hat, each=p1), z4.hat[id.k])
        
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
        pl.1000<-rbind(pl.1000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.1000.sum<-colSums(pl.1000)/n
    
  } else{ pl.1000.sum <- 0 }
  
  
  pl.0000 <- matrix(rep(0,p),1)
  rm.1i <- which(ord.res[idx.0000,1]==n1)
  rm.2i <- which(ord.res[idx.0000,2]==n2)
  rm.3i <- which(ord.res[idx.0000,3]==n3)
  rm.4i <- which(ord.res[idx.0000,4]==n4)
  if( (length(rm.1i)+length(rm.2i)+length(rm.3i)+length(rm.4i))>0) idx.0000 <- idx.0000[-unique(c(rm.1i,rm.2i,rm.3i,rm.4i))]
  if(length(idx.0000)>0L){
    for (j in idx.0000){
      idx.h <- ord.res[j, 1]
      idx.v <- ord.res[j, 2]
      idx.l <- ord.res[j, 3]
      idx.k <- ord.res[j, 4]
      
      z1.hat <- est.z1.r[(idx.h+1):n1, j]
      z2.hat <- est.z2.r[(idx.v+1):n2, j]
      z3.hat <- est.z3.r[(idx.l+1):n3, j]
      z4.hat <- est.z4.r[(idx.k+1):n4, j]
      
      p1 <- length(z1.hat)
      p2 <- length(z2.hat)
      p3 <- length(z3.hat)
      p4 <- length(z4.hat)
      
      nt3 <- (idx.l+1):n3
      nt4 <- (idx.k+1):n4
      
      nfyz123 <- matrix(rep(0,p),1)
      dfyz123 <- 0
      
      for( id.kk in 1:p4){
        
        kk <- nt4[id.kk]
        
        for( id.k in 1:p3){
          
          k <- nt3[id.k]
          prob<-f12[(idx.h+1):n1, (idx.v+1):n2]*f3[k]*f4[kk]
          
          x.s <- outer(rep(1,p1*p2), X[j,])
          XZ123 <- cbind(x.s, rep(z1.hat, p2), rep(z2.hat, each=p1), z3.hat[id.k], z4.hat[id.kk])
          
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
      }
      
      if(dfyz123!=0){
        pl.0000<-rbind(pl.0000, nfyz123/dfyz123)
      }
      x.s<-NULL
      XZ123<-NULL
      prob<-NULL
      Prob<-err.yz123<-fyz123<-NULL
      z1.hat<-NULL
      z2.hat<-NULL
      z3.hat<-NULL
      z4.hat<-NULL
    }
    
    pl.0000.sum<-colSums(pl.0000)/n
    
  } else{ pl.0000.sum <- 0 }
  
  pl.1111.sum+pl.0111.sum+pl.1011.sum+pl.1101.sum+pl.1110.sum+
    pl.0011.sum+pl.0101.sum+pl.0110.sum+pl.1001.sum+pl.1010.sum+pl.1100.sum+
    pl.0001.sum+pl.0010.sum+pl.0100.sum+pl.1000.sum+pl.0000.sum
  
}

