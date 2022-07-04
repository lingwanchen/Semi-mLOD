#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

double dnorm_cust(double x, double m, double s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

struct nfyz1_calculator : public Worker {
  const double err_yz1_0;
  const int qq;
  const RVector<int> ind_k;
  const RVector<double> theta;
  const RMatrix<double> zk;
  const RMatrix<double> pmat;
  const RVector<int> pdim;
  const double sigma0;
  const int np;
  const std::vector<long unsigned int> dim_prod;

  std::vector<double> out_nfyz1;

  nfyz1_calculator(const double err_yz1_0, const int qq,
		   const IntegerVector ind_k, const NumericVector theta,
		   const NumericMatrix zk, const NumericMatrix pmat,
		   const IntegerVector pdim, const double sigma0,
		   const int np,
		   const std::vector<long unsigned int> dim_prod) :
    		   err_yz1_0(err_yz1_0), qq(qq), ind_k(ind_k), theta(theta), zk(zk), pmat(pmat), 
			   pdim(pdim), sigma0(sigma0), np(np), dim_prod(dim_prod), out_nfyz1() {
				     out_nfyz1.resize(np+2, 0);
				   }
  nfyz1_calculator(const nfyz1_calculator& nfyz1_calculator, Split) :
    err_yz1_0(nfyz1_calculator.err_yz1_0), qq(nfyz1_calculator.qq),
    ind_k(nfyz1_calculator.ind_k), theta(nfyz1_calculator.theta),
    zk(nfyz1_calculator.zk), pmat(nfyz1_calculator.pmat),
    pdim(nfyz1_calculator.pdim), sigma0(nfyz1_calculator.sigma0),
    np(nfyz1_calculator.np), dim_prod(nfyz1_calculator.dim_prod),
    out_nfyz1() {
		out_nfyz1.resize(np+2, 0);
  }

  void operator()(std::size_t begin, std::size_t end) {
    double err_yz1;
    double dn_err_yz1;
    std::size_t i;
    std::size_t cur_i = begin;
    int j;
    int cur_coords[np];
    double prob = 1;

    for (j=np-1; j>0; j--) {
      cur_coords[j] = cur_i / dim_prod[j-1];
      cur_i %= dim_prod[j-1];
      prob *= pmat(cur_coords[j], j);
    }
	
    cur_coords[0] = cur_i;
    prob *= pmat(cur_coords[0], 0);
	
    for (i=begin; i<end; i++) {
      err_yz1 = err_yz1_0;
	  
      for (j=0; j<np; j++) {
		  err_yz1 -= theta[qq+ind_k[j]-1] * zk(j, cur_coords[j]);
      }
	  
      dn_err_yz1 = dnorm_cust(err_yz1, 0, sigma0);
      out_nfyz1[0] += dn_err_yz1 * prob * err_yz1;
	  
      for (j=1; j<=np; j++) {
		  out_nfyz1[j] += dn_err_yz1 * prob * err_yz1 * zk(j-1, cur_coords[j-1]);
      }
	  
      out_nfyz1[np+1] += dn_err_yz1 * prob;
      cur_coords[0]++;
      j = 0;
	  
      while (cur_coords[j]==pdim[j] && j<np) {
		  cur_coords[j] = 0;
		  prob *= pmat(0, j) / pmat(pdim[j]-1, j);
		  j++;
		  cur_coords[j]++;
      }
	  
      if (j<np) {
		  prob *= pmat(cur_coords[j],j) / pmat(cur_coords[j]-1,j);
      }
    }
  }

  void join(const nfyz1_calculator& rhs) {
    int i;

    for (i=0; i<(np+2); i++) {
      out_nfyz1[i] += rhs.out_nfyz1[i];
    }
  }
};

struct nfyz1_calculator_mc : public Worker {
  const double err_yz1_0;
  const int qq;
  const RVector<int> ind_k;
  const RVector<double> theta;
  const RMatrix<double> zk;
  const RMatrix<double> pmat;
  const double sigma0;
  const int np;
  const RMatrix<int> sampmat;

  std::vector<double> out_nfyz1;

  nfyz1_calculator_mc(const double err_yz1_0, const int qq,
		      const IntegerVector ind_k, const NumericVector theta,
		      const NumericMatrix zk, const NumericMatrix pmat,
		      const double sigma0, const int np,
		      const IntegerMatrix sampmat) :
    			  err_yz1_0(err_yz1_0), qq(qq), ind_k(ind_k), theta(theta), zk(zk), pmat(pmat), 
				  sigma0(sigma0), np(np), sampmat(sampmat), out_nfyz1() {
					  out_nfyz1.resize(np+2, 0);
					  }
  nfyz1_calculator_mc(const nfyz1_calculator_mc& nfyz1_calculator_mc, Split) :
    err_yz1_0(nfyz1_calculator_mc.err_yz1_0), qq(nfyz1_calculator_mc.qq),
    ind_k(nfyz1_calculator_mc.ind_k), theta(nfyz1_calculator_mc.theta),
    zk(nfyz1_calculator_mc.zk), pmat(nfyz1_calculator_mc.pmat),
    sigma0(nfyz1_calculator_mc.sigma0), np(nfyz1_calculator_mc.np),
    sampmat(nfyz1_calculator_mc.sampmat), out_nfyz1() {
		out_nfyz1.resize(np+2, 0);
  }

  void operator()(std::size_t begin, std::size_t end) {
    double err_yz1;
    double dn_err_yz1;
    std::size_t i;
    int j;
    int cur_coords[np];
    double prob;

    for (i=begin; i<end; i++) {
      prob = 1;
      for (j=0; j<np; j++) {
		  cur_coords[j] = sampmat(j,i);
		  prob *= pmat(cur_coords[j], j);
      }
	  
      err_yz1 = err_yz1_0;
      for (j=0; j<np; j++) {
		  err_yz1 -= theta[qq+ind_k[j]-1] * zk(j, cur_coords[j]);
      }
	  
      dn_err_yz1 = dnorm_cust(err_yz1, 0, sigma0);
      out_nfyz1[0] += dn_err_yz1 * prob * err_yz1;
	  
      for (j=1; j<=np; j++) {
		  out_nfyz1[j] += dn_err_yz1 * prob * err_yz1 * zk(j-1, cur_coords[j-1]);
      }
      out_nfyz1[np+1] += dn_err_yz1 * prob;
    }
  }

  void join(const nfyz1_calculator_mc& rhs) {
    int i;

    for (i=0; i<(np+2); i++) {
      out_nfyz1[i] += rhs.out_nfyz1[i];
    }
  }
};

// [[Rcpp::export]]
NumericMatrix calc_zk(int cenp, int ss, IntegerVector subn,
		      IntegerVector ind_k, IntegerMatrix ord_res,
		      List est_z_r, IntegerVector pdim) {
  IntegerVector dm(cenp);
  int i;
  int j;
  int max_drop;
  NumericMatrix outMA(cenp, max(pdim));

  for (i=0; i<cenp; i++) {
    for (j=0; j<cenp; j++) {
      if (j<i) {
		  dm[j] = j+2;
      }
      else if (j==i) {
		  dm[j] = 1;
      }
      else if (j>i) {
		  dm[j] = j+1;
      }
    }
    NumericMatrix temp1 = est_z_r[ind_k[i]-1];
    NumericVector temp2 = temp1(_, subn[ss-1]-1);
    max_drop = ord_res(subn[ss-1]-1, ind_k[i]-1);
    if (max_drop!=0) {
      NumericVector temp2a(temp2.size()-max_drop);
      for (j=0; j<temp2a.size(); j++) {
		  outMA(i,j) = temp2[j+max_drop];
      }
    }
    else {
      for (j=0; j<temp2.size(); j++) {
		  outMA(i,j) = temp2[j];
      }
    }
  }
  return(outMA);
}

// [[Rcpp::export]]
NumericVector calc_nfyz1(double err_yz1_0, int qq, IntegerVector ind_k,
			 NumericVector theta, NumericMatrix zk,
			 NumericMatrix pmat, IntegerVector pdim,
			 double sigma0) {
  int np = pdim.size();
  std::vector<long unsigned int> dim_prod;
  int i;

  dim_prod.resize(np, 1);
  dim_prod[0] = pdim[0];
  for (i=1; i<np; i++) {
    dim_prod[i] = pdim[i] * dim_prod[i-1];
  }

  nfyz1_calculator cur_nfyz1_calc(err_yz1_0, qq, ind_k, theta, zk, pmat, pdim, sigma0, np, dim_prod);

  parallelReduce(0, dim_prod[np-1], cur_nfyz1_calc);

  return(wrap(cur_nfyz1_calc.out_nfyz1));
}

// [[Rcpp::export]]
NumericVector calc_nfyz1_mc(double err_yz1_0, int qq, IntegerVector ind_k,
			    NumericVector theta, NumericMatrix zk,
			    NumericMatrix pmat, double sigma0,
			    IntegerVector pdim, int nmc) {
  int np = pmat.ncol();
  IntegerMatrix sampmat(np, nmc);
  int i;
  int j;
  NumericVector unifv = Rcpp::runif(np*nmc);

  for (i=0; i<np; i++) {
    for (j=0; j<nmc; j++) {
      sampmat(i,j) = (int) floor(pdim(i)*unifv[i*nmc+j]);
    }
  }

  nfyz1_calculator_mc cur_nfyz1_calc(err_yz1_0, qq, ind_k, theta, zk, pmat, sigma0, np, sampmat);

  parallelReduce(0, nmc, cur_nfyz1_calc);

  return(wrap(cur_nfyz1_calc.out_nfyz1));
}
