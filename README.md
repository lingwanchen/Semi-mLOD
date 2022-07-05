# Semi-mLOD
Proposed methods in "Semiparametric analysis of a generalized linear model with multiple covariates subject to detection limits". The implemented methods work for the continuous, binary or count type of outcome. The rank-based method and the least-square method for semiparametric AFT models are available to use. The joint distribution of residuals are estimated by either the mulitivariate K-M estimator or the marginal approximation approach. When the number of dimension is greater 3 and outcome type is continuous, the marginal approximation approach, incorprating the MCMC alogrithm in the estimation of theta, is avaiable to use.

reqfuns.R: required functions for the implementation 

aft_lod.cpp: required function for the MCMC

main-funs-mKM.R: main function for the proposed methods, including multivaritet KM estimator and marginal approximation approach (p<4). 

main-funs-LM-mKM.R: main function for the proposed method using marginal approximation approach and MCMC alogrithm where type of outcome is continuous (p>3).

R-examples.R: five examples are provided to illustrate the proposed methods along with the sample datasets (p=2, 3, 4)

Note: User also needs seven required packages ("rootSolve", "quantreg" and "SparseM", "survival", "aftgee", "Rcpp", "RcppParallel) from R before using our method
