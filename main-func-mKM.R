main.func.mKM<-function(Y, Ytype, Y.cov, Z, lod, aft.cov, trans, alpha.method, eta.method){
  # Y: outcome of interest
  # Ytype ("cts", "bin" or "count"): type of outcome Y
  # X: fully observed covariates (input should be a model-design matrix)
  # Z: is a nxp matrix for covariates subject to LOD (must be continuous)
  # lod: LOD values for Z
  # aft.cov: covariates for AFT model adjustment
  # transforz ("ne" or "exp-ne"): transform function h(t)=z: -t "ne" or exp(-t) "exp-ne"
  # alpha.method ("rank" or "LS"): rank-based with Gehan weight "rank" or least-square method with exchangeable correlation "LS"
  # eta.method ("mKM" or"marg"): multivariate-KM "mKM" or marginal "marg"
  ##  NOTE: this function works for p<4
  
   
  pp<-dim(Z)[2]
  
  if(pp>3){ output<-"not allow p>3" }else{

      if( Ytype=="cts"|| Ytype=="bin" || Ytype=="count"){
        if( alpha.method=="rank"|| alpha.method=="LS"){
          if( eta.method=="mKM"|| eta.method=="marg"){
            
        n<-length(Y)
        
        data<-as.data.frame(cbind(Y,Y.cov))
        data$X<-Y.cov  
        data$Y<-Y
        
        varname<-colnames(data$X)
        aft.var<-colnames(aft.cov)
        data$XX <- as.matrix(cbind(1, data$X))
        
        
        bbb<-varname[1]
        for( j in 2:length(varname)) bbb<-paste0(bbb,"+", varname[j])
        form.y.cov <- paste0("Y ~ ",bbb)
        
        bbb2<-aft.var[1]
        for( j in 2:length(aft.var)) bbb2<-paste0(bbb2,":margin+", aft.var[j])
        form.aft.cov <- paste0("Surv(expZ, Delta) ~ ",bbb2,":margin + margin - 1")
        
        if(trans=="ne") {gf <- nega; inv.gf <-inv.nega;}
        if(trans=="exp-ne") {gf <- exp.ne; inv.gf <-inv.exp.ne;}
        
        
        if(pp==1){
          
          data$z1<-Z[,1]
          
          l1<-lod[1]
          
          data$delta1 <- ifelse(data$z1>=l1, 1, 0)
          
          data$z1st <- pmax(data$z1, l1)
          
          data$Z <- data$z1st
          
          data$XZ <- as.matrix(cbind(data$XX, data$Z)) 
          
          idx.1 <<- which(data$delta1==1)
          idx.0 <<- which(data$delta1==0) 
          
          T1 <- inv.gf(data$z1st)
          alpha1 <- coef.rankAFT(y=T1, x=data.frame(aft.cov),  delta = data$delta1)
          r1 <- c(T1-alpha1%*%t(cbind(1,aft.cov)))
          S1.all <- survf(r1,data$delta1)
          ff <<- S1.all$prob
          n1 <<- length(ff)
          res1<- sort(unique(r1))
          
          tmp1 <- matrix(rep(r1, each = n1), n1, n)
          ttt <- ifelse(tmp1==res1, 1, 0)
          ord.res<<- n1-colSums(apply(ttt,2,cumsum))+1
          
          ttt <- tmp1 <- NULL
          
          est.z1.r <<- gf(outer(res1, c(alpha1%*%t(as.matrix(cbind(1,aft.cov)))),"+"))          
          
          form.z <-"+z1st"
          formula.cc <- as.formula(paste0(form.y.cov, form.z))
          
          if(Ytype=="cts"){
            CC.reg <- lm(formula.cc, data = data[-idx.0, ])
            sigma.0 <- summary(CC.reg)$sigma
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.LM.p1, theta, parms = list(data=data, sigma.0=sigma.0, p=p))$root)
          }
          if(Ytype=="bin"){
            CC.reg <- glm(formula.cc, family =binomial(link="logit"), data = data[-idx.0, ])
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.BN.p1, theta, parms = list(data=data, p=p))$root)
          }
          if(Ytype=="count"){
            CC.reg <- glm(formula.cc,  family=poisson(link="log"), data = data[-idx.0, ])
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.PO.p1, theta, parms = list(data=data, p=p))$root)
          }
          
          rm("idx.1","idx.0", pos = ".GlobalEnv")
          rm("ff", "n1",pos = ".GlobalEnv")
          rm("ord.res","est.z1.r", pos = ".GlobalEnv")
          
        }
        
        if(pp==2){
          
          data$z1<-Z[,1]
          data$z2<-Z[,2]
          
          l1<-lod[1]
          l2<-lod[2]
          
          data$delta1 <- ifelse(data$z1>=l1, 1, 0)
          data$delta2 <- ifelse(data$z2>=l2, 1, 0)
          
          data$z1st <- pmax(data$z1, l1)
          data$z2st <- pmax(data$z2, l2)
          
          data$Z <- cbind(data$z1st, data$z2st)
          
          data$XZ <- as.matrix(cbind(data$XX, data$Z))  
          idx.c1 <- which(data$delta1==0)
          idx.c2 <- which(data$delta2==0)
          idx.cen <- unique(c(idx.c1, idx.c2))  
          
          idx.11 <<- which(data$delta1*data$delta2==1)
          idx.10 <<- which(data$delta1*(1-data$delta2)==1)
          idx.01 <<- which((1-data$delta1)*data$delta2==1)
          idx.00 <<- which((1-data$delta1)*(1-data$delta2)==1)
          
          
          
          T1 <- inv.gf(data$z1st)
          T2 <- inv.gf(data$z2st)
          
          if(alpha.method=="rank"){
            alpha1 <- coef.rankAFT(y=T1, x=data.frame(aft.cov),  delta = data$delta1)
            alpha2 <- coef.rankAFT(y=T2, x=data.frame(aft.cov),  delta = data$delta2)
          }
          
          if(alpha.method=="LS"){
            data$id<-1:n
            alpha.x<-data.frame(t1=exp(T1), t2=exp(T2), delta1=data$delta1, delta2=data$delta2, aft.cov, id=data$id)
            long.my<-reshape(alpha.x, idvar="id", direction="long", varying=list(t=c(1,2), delta=c(3,4)), v.names = c("expZ", "Delta") )
            long.my$margin<-factor(long.my$time)
            
            fit1.ex <- aftgee(as.formula(form.aft.cov), id = id, margin=margin, data = long.my, corstr = "ex", B=0)$coefficients
            alpha1 <- fit1.ex[c(1,3,5),2]
            alpha2 <- fit1.ex[c(2,4,6),2]
            alphax.x <- long.my <- fit1.ex <- NULL
          }
          
          r1 <- c(T1-alpha1%*%t(cbind(1,aft.cov)))
          r2 <- c(T2-alpha2%*%t(cbind(1,aft.cov)))
          
          if(eta.method=="mKM"){
            DD <- cbind(r1, r2, data$delta1, data$delta2)
            BS <- DB.BS(DD)
            SS<-BS$BS
            S12<-outer(BS$S1,BS$S2,"+")
            FF<-SS+1-S12
            FF <- apply(FF,2,cummax)
            FF <- t(apply(FF,1,cummax))
            n1 <<- dim(FF)[1]
            n2 <<- dim(FF)[2]
            
            res1<-BS$t1
            res2<-BS$t2
            ff <<- FF + cbind(0,rbind(0,FF[1:(n1-1),1:(n2-1)])) - rbind(0,FF[1:(n1-1),]) - cbind(0,FF[,1:(n2-1)])
          }
          if(eta.method=="marg"){
            S1.all <- survf(r1,data$delta1)
            S2.all <- survf(r2,data$delta2)
            res1<-S1.all$time[S1.all$prob>0]
            res2<-S2.all$time[S2.all$prob>0]
            n1 <<- length(res1)
            n2 <<- length(res2)
            
            f1F2 <<- matrix(c(S1.all$prob[S1.all$prob>0]),n1,n2,byrow=F)
            f2F1 <<- matrix(c(S2.all$prob[S2.all$prob>0]),n2,n1,byrow=F)
            ff <<- f1F2[,1]%*%t(f2F1[,1])
            
          }
          
          tmp1 <- matrix(rep(r1, each = n1), n1, n)
          ord1 <- colSums(ifelse(tmp1>=res1, 1, 0))
          tmp2 <- matrix(rep(r2, each = n2), n2, n)
          ord2 <- colSums(ifelse(tmp2>=res2, 1, 0))
          ttt <- tmp1 <- tmp2 <- NULL
          ord.res <<-cbind(ord1,ord2)
          est.z1.r <<- gf(outer(res1, c(alpha1%*%t(as.matrix(cbind(1,aft.cov)))),"+"))
          est.z2.r <<- gf(outer(res2, c(alpha2%*%t(as.matrix(cbind(1,aft.cov)))),"+"))
          
          form.z <-"+z1st+z2st"
          formula.cc <- as.formula(paste0(form.y.cov, form.z))
          
          if(Ytype=="cts"){
            CC.reg <- lm(formula.cc, data = data[-idx.cen, ])
            sigma.0 <- summary(CC.reg)$sigma
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.LM.p2, theta, parms = list(data=data, sigma.0=sigma.0, p=p, n=n))$root)
          }
          if(Ytype=="bin"){
            CC.reg <-  glm(formula.cc, family =binomial(link="logit"), data = data[-idx.cen, ])
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.BN.p2, theta, parms = list(data=data, p=p, n=n))$root)
            
          }
          if(Ytype=="count"){
            CC.reg <-  glm(formula.cc,  family=poisson(link="log"), data = data[-idx.cen, ])
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.PO.p2, theta, parms = list(data=data, p=p, n=n))$root)
            
          }
          rm("idx.11","idx.10","idx.01","idx.00" , pos = ".GlobalEnv")
          rm("ff", "n1","n2", pos = ".GlobalEnv")
          rm("ord.res","est.z1.r","est.z2.r", pos = ".GlobalEnv")
          
        }
        
        if(pp==3){
          
          data$z1<-Z[,1]
          data$z2<-Z[,2]
          data$z3<-Z[,3]
          
          l1<-lod[1]
          l2<-lod[2]
          l3<-lod[3]
          
          data$delta1 <- ifelse(data$z1>=l1, 1, 0)
          data$delta2 <- ifelse(data$z2>=l2, 1, 0)
          data$delta3 <- ifelse(data$z3>=l3, 1, 0)
          
          data$z1st <- pmax(data$z1, l1)
          data$z2st <- pmax(data$z2, l2)
          data$z3st <- pmax(data$z3, l3)
          
          data$Z <- cbind(data$z1st, data$z2st, data$z3st)
          
          data$XZ <- as.matrix(cbind(data$XX, data$Z)) 
          
          idx.c1 <- which(data$delta1==0)
          idx.c2 <- which(data$delta2==0)
          idx.c3 <- which(data$delta3==0)
          idx.cen <- unique(c(idx.c1, idx.c2, idx.c3))
          
          idx.111 <<- which(data$delta1*data$delta2*data$delta3==1)
          idx.100 <<- which(data$delta1*(1-data$delta2)*(1-data$delta3)==1)
          idx.010 <<- which(data$delta2*(1-data$delta1)*(1-data$delta3)==1)
          idx.001 <<- which(data$delta3*(1-data$delta1)*(1-data$delta2)==1)
          idx.110 <<- which(data$delta1*data$delta2*(1-data$delta3)==1)
          idx.101 <<- which(data$delta1*data$delta3*(1-data$delta2)==1)
          idx.011 <<- which(data$delta2*data$delta3*(1-data$delta1)==1)
          idx.000 <<- which((1-data$delta1)*(1-data$delta2)*(1-data$delta3)==1)
          
          T1 <- inv.gf(data$z1st)
          T2 <- inv.gf(data$z2st)
          T3 <- inv.gf(data$z3st)    
          
          if(alpha.method=="rank"){
            alpha1 <- coef.rankAFT(y=T1, x=data.frame(aft.cov),  delta = data$delta1)
            alpha2 <- coef.rankAFT(y=T2, x=data.frame(aft.cov),  delta = data$delta2)
            alpha3 <- coef.rankAFT(y=T3, x=data.frame(aft.cov),  delta = data$delta3)
          }
          
          if(alpha.method=="LS"){
            data$id<-1:n
            alpha.x<-data.frame(t1=exp(T1), t2=exp(T2), t3=exp(T3), delta1=data$delta1, delta2=data$delta2, delta3=data$delta3, aft.cov, id=data$id)
            long.my<-reshape(alpha.x, idvar="id", direction="long", varying=list(t=c(1,2,3), delta=c(4,5,6)), v.names = c("expZ", "Delta") )
            long.my$margin<-factor(long.my$time)
            
            fit1.ex <- aftgee(as.formula(form.aft.cov), id = id, margin=margin, data = long.my, corstr = "ex", B=0)$coefficients
            alpha1 <- fit1.ex[c(1,4,7),2]
            alpha2 <- fit1.ex[c(2,5,8),2]
            alpha3 <- fit1.ex[c(3,6,9),2]
            alphax.x <- long.my <- fit1.ex <- NULL
          }
          
          r1 <- c(T1-alpha1%*%t(cbind(1,aft.cov)))
          r2 <- c(T2-alpha2%*%t(cbind(1,aft.cov)))
          r3 <- c(T3-alpha3%*%t(cbind(1,aft.cov)))
          
          if(eta.method=="mKM"){
            DD <- cbind(r1, r2, r3, data$delta1, data$delta2, data$delta3)
            TSF <- TriKM(DD)
            
            S123<-TSF$TSF
            dp <- dim(S123)
            
            S1<-TSF$S1
            S2<-TSF$S2
            S3<-TSF$S3
            S12<-TSF$S12
            S13<-TSF$S13
            S23<-TSF$S23
            res1 <- TSF$T1 
            res2 <- TSF$T2 
            res3 <- TSF$T3
            
            TSF<-NULL
            
            One.m <- array(1, dp)
            S12.m <- array(S12, dp)
            S13.m <- aperm(array(S13, dp[c(1,3,2)]), c(1,3,2))
            S23.m <- aperm(array(S23, dp[c(2,3,1)]), c(3,1,2))
            S1.m <- array(rep(S1,dp[2]*dp[3]), dp)
            S2.m <- aperm(array(rep(S2,dp[1]*dp[3]), dp[c(2,1,3)]), c(2,1,3))
            S3.m <- aperm(array(rep(S3,dp[1]*dp[2]), dp[c(3,1,2)]), c(2,3,1))
            
            FFF <- One.m+S12.m+S13.m+S23.m-S1.m-S2.m-S3.m-S123
            
            one.m<-S12.m<-S13.m<-S23.m<-S1.m<-S2.m<-S3.m<-S123<-NULL
            
            FFF <- apply(FFF,c(2,3),cummax)
            FFF <- aperm(apply(FFF,c(1,3),cummax),c(2,1,3))
            FFF <- aperm(apply(FFF,c(1,2),cummax),c(2,3,1))
            
            n1 <<- length(res1)
            n2 <<- length(res2)
            n3 <<- length(res3)

            f123<<- FFF[-1,-1,-1] - FFF[-1,-1,-dp[3]] - FFF[-1,-dp[2],-1] - FFF[-dp[1],-1,-1]+
              FFF[-1,-dp[2],-dp[3]] + FFF[-dp[1],-1,-dp[3]] + FFF[-dp[1],-dp[2],-1] - FFF[-dp[1],-dp[2],-dp[3]]
            
          }
          if(eta.method=="marg"){
            S1.all <- survf(r1,data$delta1)
            S2.all <- survf(r2,data$delta2)
            S3.all <- survf(r3,data$delta3)
            
            res1<-S1.all$time[S1.all$prob>0]
            res2<-S2.all$time[S2.all$prob>0]
            res3<-S3.all$time[S3.all$prob>0]  
            
            n1 <<- length(res1)
            n2 <<- length(res2)
            n3 <<- length(res3)
            
            f12F3 <<- array(outer(c(S1.all$prob[S1.all$prob>0]),c(S2.all$prob[S2.all$prob>0]),"*"),c(n1,n2,n3))
 
            f123<<-outer(f12F3[,,1],c(S3.all$prob[S3.all$prob>0]),"*")
            
            S3.all<-S1.all<-S2.all<-NULL
            
          }
          
          FFF<-NULL
          
          tmp1 <- matrix(rep(r1, each = n1), n1, n) 
          ord1 <- colSums(ifelse(tmp1>=res1, 1, 0))
          
          tmp2 <- matrix(rep(r2, each = n2), n2, n) 
          ord2 <- colSums(ifelse(tmp2>=res2, 1, 0))
          
          tmp3 <- matrix(rep(r3, each = n3), n3, n) 
          ord3 <- colSums(ifelse(tmp3>=res3, 1, 0))
          
          ttt <- tmp1 <- tmp2 <-tmp3 <- NULL
          
          ord.res <<-cbind(ord1,ord2,ord3)
          
          est.z1.r <<- gf(outer(res1, c(alpha1%*%t(as.matrix(cbind(1,aft.cov)))),"+"))
          est.z2.r <<- gf(outer(res2, c(alpha2%*%t(as.matrix(cbind(1,aft.cov)))),"+"))
          est.z3.r <<- gf(outer(res3, c(alpha3%*%t(as.matrix(cbind(1,aft.cov)))),"+"))
          
          form.z <-"+z1st+z2st+z3st"
          formula.cc <- as.formula(paste0(form.y.cov, form.z))
          
          if(Ytype=="cts"){
            CC.reg <- lm(formula.cc, data = data[-idx.cen, ])
            sigma.0 <- summary(CC.reg)$sigma
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.LM.p3, theta, parms = list(data=data, sigma.0=sigma.0, p=p, n=n))$root)
          }
          if(Ytype=="bin"){
            CC.reg <-  glm(formula.cc, family =binomial(link="logit"), data = data[-idx.cen, ])
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.BN.p3, theta, parms = list(data=data, p=p, n=n))$root)
          }
          if(Ytype=="count"){
            CC.reg <-  glm(formula.cc,  family=poisson(link="log"), data = data[-idx.cen, ])
            
            theta <- CC.reg$coefficients
            p <- length(theta)    
            
            theta.hat <- try(multiroot(pslogL.PO.p3, theta, parms = list(data=data, p=p, n=n))$root)
          }
          
          rm("idx.111","idx.110","idx.101","idx.011","idx.001","idx.010","idx.100","idx.000" , pos = ".GlobalEnv")
          rm("f123", "n1","n2","n3", pos = ".GlobalEnv")
          rm("ord.res","est.z1.r","est.z2.r","est.z3.r", pos = ".GlobalEnv") 
        }
        
        output<-theta.hat  
          }else{
            output <- "only allows multivariate K-M or marginal approach"
          }
        }else{
            output <- "only allows rank-based or least-square method"
        }
        }
      else{
        output <- "only allows cts, binary or count"
      }
  }
  
  return(list(theta.hat=output))
}
