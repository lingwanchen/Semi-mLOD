# Semi-mLOD
Proposed methods in "Semiparametric analysis of a generalized linear model with multiple covariates subject to detection limits". The implemented methods work for the continuous, binary or count type of outcome and p<5. The rank-based method and the least-square method for semiparametric AFT models are available to use. The joint distribution of residuals are estimated by either the mulitivariate K-M estimator or the marginal approach (marginal approach is available for p=4)

reqfuns.R: required functions for the implementation 

main.funs.R: main function for the proposed methods

R-examples.R: some examples are provided to illustrate the proposed methods along with the sample datasets

Note: User also needs five required packages ("rootSolve", "quantreg" and "SparseM", "survival", "aftgee") from R before using our method
