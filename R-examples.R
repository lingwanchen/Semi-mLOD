### Three examples are provided here:
### 1. p=2 and type of outcome is binary.
### 2. p=3 and type of outcome is count 
### 3. p=4 and type of outcome is continuous


source("reqfuns.r")   ## required functions
source("main-func.r") ## main function of our methods. See details in "main-func.R"

library(rootSolve)
library(quantreg)
library(SparseM)
library(survival)
library(aftgee)
library(copula)

##########################################################################################
## p=2, Ytype="bin", trans="ne", alpha.method="rank" or "LS", eta.method="mKM" or "marg" # 
##########################################################################################
data <- read.table("data-bin-p2.txt",sep=' ') 
colnames(data) <- c("Y","z1","z2","x1","x2")       

Y <- data$Y
Y.covX <- cbind(data$x1, data$x2)
colnames(Y.covX) <- c("xbin", "xcts")
Y.covZ <- cbind(data$z1, data$z2)
lod <- c(1, 1)

## alpha.method="rank" and eta.method="mKM"
(ybin.p2.rk.mkm<-main.func(Y=Y, Ytype="bin", Y.cov=Y.covX, Z=Y.covZ, lod=lod, aft.cov=Y.covX, trans="ne", alpha.method="rank", eta.method="mKM"))

## alpha.method="LS" and eta.method="marg"
(ybin.p2.ls.marg<-main.func(Y=Y, Ytype="bin", Y.cov=Y.covX, Z=Y.covZ, lod=lod, aft.cov=Y.covX, trans="ne", alpha.method="LS", eta.method="marg"))


########################################################################
## p=3, Ytype="count", trans="ne", alpha.method="LS", eta.method="mKM" #
########################################################################
data <- read.table("data-count-p3.txt",sep=' ') 
colnames(data) <- c("Y","z1","z2","z3","x1","x2")       

Y <- data$Y
Y.covX <- cbind(data$x1, data$x2)
colnames(Y.covX) <- c("xbin","xcts")
Y.covZ <- cbind(data$z1,data$z2,data$z3)
lod <- c(1.2,1.2,1.3)

(ycount.p3.mkm<-main.func(Y=Y, Ytype="count", Y.cov=Y.covX, Z=Y.covZ, lod=lod, aft.cov=Y.covX, trans="ne", alpha.method="LS", eta.method="mKM"))


#########################################################################
## p=4, Ytype="cts", trans="ne", alpha.method="rank", eta.method="marg" #
#########################################################################
data <- read.table("data-cts-p4.txt",sep=' ') 
colnames(data) <- c("Y","z1","z2","z3","z4","x1","x2","x3","x4")       

Y <- data$Y
Y.covX <- cbind(data$x1, data$x2, data$x3, data$x4)
colnames(Y.covX) <- c("xbin","xcts1","xcts2","xcts3")
Y.covZ <- cbind(data$z1,data$z2,data$z3,data$z4)
lod<-c(1.3,1.3,1.5,1.5)

(ycts.p4.rk.marg<-main.func(Y=Y, Ytype="cts", Y.cov=Y.covX, Z=Y.covZ, lod=lod, aft.cov=Y.covX, trans="ne", alpha.method="rank", eta.method="marg"))


