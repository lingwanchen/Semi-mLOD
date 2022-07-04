# This function works for the continuous outcome and p>2. For p<3, use main.func.mKM function instead.
# The marginal approach Wwith rank-based with Gehan weight for alpha estimation is used to estimate the eta (the residual term distribution).
# The MCMC algorithm was implemented in the estimation of theta.

main.func.LM.marg.mc<-function(Y, Y.cov, Z, lod, aft.cov, trans, seedmc, contl){
  # Y: outcome of interest, cts type
  # Y.cov: fully observed covariates (input should be a model-design matrix)
  # Z: is a nxp matrix for covariates subject to LOD (must be continuous)
  # lod: LOD values for Z
  # aft.cov: covariates for AFT model adjustment
  # transforz ("ne" or "exp-ne"): transform function h(t)=z: -t "ne" or exp(-t) "exp-ne"
  # seedmc: seed for MCMC
  ##  NOTE: this function works for p>2
  
  pp<-dim(Z)[2]
    
  if(pp<3){ 
    output<-"This function works for p>2" 
  }else{
    n<-length(Y)
    
    data<-as.data.frame(cbind(Y,Z,Y.cov))
    data$X<-Y.cov  
    data$Y<-Y
    
    varname<-colnames(data$X)
    aft.var<-colnames(aft.cov)
    data$XX <- as.matrix(cbind(1, data$X))
    
    
    bbb<-varname[1]
    for( j in 2:length(varname)) bbb<-paste0(bbb,"+", varname[j])
    form.y.cov <- paste0("Y ~ ",bbb)
    
    
    if(trans=="ne") {gf <- nega; inv.gf <-inv.nega;}
    if(trans=="exp-ne") {gf <- exp.ne; inv.gf <-inv.exp.ne;}
    
    cen<-data$Zst<-matrix(0, n, pp)
    
    # calculate the censoring
    # Zst is the matrix that has LOD in Z. If Z matrix includes LOD, then Zst=Z
    for( s in 1:pp){
      cen[, s]<-ifelse(data[, (s+1)]> lod[s], 1, 0)
      data$Zst[, s]<-pmax(data[, (s+1)], lod[s])
    }
    
    # indicator the subject in the particular missing group
    # miss.indx$miss0: subjects have all Z observed. ie CC
    # miss.indx$miss1: subjects with one Z below LOD, and so on...
    
    miss.indx <- c()
    miss.indx[[paste0("miss",0)]]<-which(rowSums(cen)== c(pp) )
    for( j in 1:pp){
      miss.indx[[paste0("miss",j)]]<-which(rowSums(cen)== c(pp-j) )
      cb<-combn(pp,j)
      ncb<-dim(cb)[2]
      for( jj in 1:ncb){
        aa<-paste0(cb[,jj], collapse = "")
        mm1<-miss.indx[[paste0("miss",j)]]
        if(length(mm1)>1){
          miss.indx[[paste0("p",aa)]]<-miss.indx[[paste0("miss",j)]][which( rowSums(as.matrix(cen[ miss.indx[[paste0("miss",j)]] , cb[,jj]],,j))== 0)]
        }
        if(length(mm1)==1){
          miss.indx[[paste0("p",aa)]]<-miss.indx[[paste0("miss",j)]][which( sum(as.matrix(cen[ miss.indx[[paste0("miss",j)]] , cb[,jj]],,j))== 0)]
          
        }
        
      }
    }
    
    data$XZ <- as.matrix(cbind(data$XX, data$Zst))  # for estimation
    
    ## initial theta from the regression using CC
    CC.reg <- lm(Y ~ X+Zst, data = data[ miss.indx$miss0, ])
    sigma.0 <- summary(CC.reg)$sigma
    
    theta <- CC.reg$coefficients
    p <- length(theta)
    
    num.na<-sum(is.na(theta))
    if( num.na==0 ){
      # estimate alpha in the AFT model
      alpha.rk <- TT <- resi <-c()
      j<-1
      while( j < (pp+1) ){
        TT <-cbind(TT, inv.gf(data$Zst[, j]))
        saft<-try(coef.rankAFT(y=TT[, j], x=data.frame(aft.cov),  delta =cen[, j]))
        if( class(saft)!="try-error"){
          alpha.rk <- rbind(alpha.rk, saft)
          resi<-cbind(resi, c(TT[, j]- alpha.rk[j,]%*%t(cbind(1,aft.cov))))
          j<-j+1
        }else{
          j<-pp+10
        }
      }
      
      if( class(saft)!="try-error"){ # All AFT models work
        
        ## estimate the KM estimator for each residual term (T-hat.T)
        ## for each subj, estimate hat.Z(s), eta(s), the corresponding indicator and the lower bound in the integration
        S.time <- c()
        res.or <- c()
        prob.or<- c()
        n.v <- c()
        ord.res <-c()
        est.z.r<-c()
        ind.use<-c()
        for( j in 1:pp){
          S.time[[paste0("S",j)]]<-survf(resi[, j], cen[, j])
          res.or[[paste0("res",j)]]<-S.time[[paste0("S",j)]]$time[S.time[[paste0("S",j)]]$prob>0]
          prob.or[[paste0("p",j)]]<-S.time[[paste0("S",j)]]$prob[S.time[[paste0("S",j)]]$prob>0]
          nk<-length(prob.or[[paste0("p",j)]])
          
          n.v<-c(n.v, nk)
          
          tmp <- matrix(rep(resi[, j], each = nk), nk, n)
          ord1 <- colSums(ifelse(tmp>=res.or[[paste0("res",j)]], 1, 0))
          
          ind.use[[paste0("Z",j)]] <-ifelse(tmp>=res.or[[paste0("res",j)]], 0, 1)
          
          ord.res<-cbind(ord.res, ord1)
          
          est.z.r[[paste0("Z",j)]] <- gf(outer(res.or[[paste0("res",j)]], c(alpha.rk[j,]%*%t(as.matrix(cbind(1,aft.cov)))),"+"))*ind.use[[paste0("Z",j)]]
        }
        
        ## the parameters required in the optima.
        parms = list(data=data, miss.indx=miss.indx, est.z.r=est.z.r, ind.use=ind.use, res.or=res.or, 
                     prob.or=prob.or, n.v=n.v, ord.res=ord.res, sigma.0=sigma.0, seedmc=seedmc)
        
        (theta.hat <- try(multiroot(pslogL.LM.marg.mc, theta, parms = parms, ctol=contl)$root))
        
       # rm("est.z.r","ind.use","res.or","prob.or","n.v","ord.res","tmp","S.time", pos = ".GlobalEnv")

        
        output<-theta.hat
        
      }else{
        TT<-inv.gf(data$Zst[, 1])
        theta.hat<-try(coef.rankAFT(y=TT[1:2], x=data.frame(data$X[1:2,]),  delta =cen[1:2, 1]))
        output<-theta.hat
      }
      
    }else{
      TT<-inv.gf(data$Zst[, 1])
      theta.hat<-try(coef.rankAFT(y=TT[1:2], x=data.frame(data$X[1:2,]),  delta =cen[1:2, 1]))
      output<-theta.hat
    }
    
   # rm("miss.indx","cen", pos = ".GlobalEnv")
    
  }
  
  
  return(list(theta.hat=output))
}
