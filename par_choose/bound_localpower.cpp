// This file is a small alteration of the bound.cpp file.  For the confidence 
// band experiment, it produces sup-norm and L2-norm statistics of only the 
// lower bound.
// 
// This file contains functions used to conduct simulation experiements 
// and generally to conduct inference for first and second order PT 
// stochastic dominance, whether distributions are point- or 
// partially-identified.

#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// Collect a pooled sample for the difference empirical process, used
// for FOSD process.
// [[Rcpp::export]]
NumericVector deval(NumericVector sam1, NumericVector sam0) {
  std::vector<double> eval;
  eval.insert(eval.end(), sam1.begin(), sam1.end());
  eval.insert(eval.end(), sam0.begin(), sam0.end());
  std::sort(eval.begin(), eval.end());
  std::vector<double>::iterator it;
  it = std::unique(eval.begin(), eval.end());
  eval.resize(std::distance(eval.begin(), it));
  return wrap(eval);
}

// An ecdf() function
// eval is where to evaluate the empirical CDF 
// sobs is sorted observations.
// [[Rcpp::export]]
NumericVector edf(NumericVector eval, NumericVector sobs) {
  int neval = eval.size();
  IntegerVector ans(neval);
  for (int i = 0; i < neval; ++i)
    ans[i] = std::upper_bound(sobs.begin(), sobs.end(), eval[i]) - sobs.begin();
  return as<NumericVector>(ans) / sobs.size();
}

// Compute the L2 norm of a one-dimensional process.
// fun is the function, evaluated at eval
// [[Rcpp::export]]
double l2norm(NumericVector eval, NumericVector fun) {
  int pm1 = eval.size() - 1;
  NumericVector inside = no_init(pm1);
  for (int i = 0; i < pm1; ++i) {
    inside(i) = 0.5 * (pow(fun(i), 2.0) + pow(fun(i+1), 2.0));
    inside(i) *= (eval(i + 1) - eval(i));
  }
  return sqrt(sum(inside));
}

// Compute lower and upper bounds for the CDF of the difference 
// between two samples.
// eval is where to evaluate it (some sort of grid usually).
// First row of bounds is lower bound, second row is upper bound.

// [[Rcpp::export]]
NumericMatrix lohi(NumericVector eval, NumericVector sam1, NumericVector sam0) {
  int ne = eval.size();
  std::sort(sam1.begin(), sam1.end());
  std::sort(sam0.begin(), sam0.end());
  NumericVector F1 = edf(eval, sam1);
  NumericMatrix bounds = no_init(2, ne);
  NumericVector shift0(sam0.size());
  NumericVector F0(ne);
  NumericVector dpro(ne);
  for (int i = 0; i < ne; ++i) {
    for (int j = 0; j < sam0.size(); ++j)
      shift0[j] = sam0[j] + eval[i];
    F0 = edf(eval, shift0);
    dpro = F1 - F0;
    dpro.push_back(0.0);
    bounds(0, i) = max(dpro);
    bounds(1, i) = 1.0 + min(dpro);
  }
  return bounds;
}

// This function estimates some sets that are used in the bootstrap 
// algorithm later to make directional derivative estimates.  The sets
// are epsilon-maximizers and minimizers (sets marked with 1).

// [[Rcpp::export]]
List epsmax_cb(NumericVector eval, NumericVector sam1, NumericVector sam0, 
              double eps) {
  int ne = eval.size();
  std::sort(sam1.begin(), sam1.end());
  std::sort(sam0.begin(), sam0.end());
  NumericVector F1 = edf(eval, sam1);
  NumericMatrix extreme = no_init(2, ne);
  NumericVector phiL(ne);
  NumericVector phiU(ne);
  NumericVector shift0(sam0.size());
  NumericVector F0(ne);
  NumericVector dpro(ne);
  IntegerMatrix epsmax(ne, ne);
  IntegerMatrix epsmin(ne, ne);
  for (int i = 0; i < ne; ++i) {
    for (int j = 0; j < sam0.size(); ++j)
      shift0[j] = sam0[j] + eval[i];
    F0 = edf(eval, shift0);
    dpro = F1 - F0;
    extreme(0, i) = max(dpro);
    extreme(1, i) = min(dpro);
    for (int j = 0; j < ne; ++j) {
      if (dpro[j] >= extreme(0, i) - eps) epsmax(j, i) = 1;
      if (dpro[j] <= extreme(1, i) + eps) epsmin(j, i) = 1;
    }
  }
  return List::create(_["epsmaxL"] = epsmax, _["epsminU"] = epsmin);
}

// Uniform confidence bands for the functions that bound the CDF of the 
// treatment effect when the dependency between pre- and post-treatment
// variables is unknown.  eval is where to evaluate the functions, sam0 and sam1
// are samples of observations, R is the number of bootstrap repetitions and eps
// is a small number used to estimate epsilon-maximizers in one argument for
// each value of the other argument.  The function returns two bootstrap
// samples, one for each of the upper and lower bound functions, where the
// observations are all supremum norm statistics.

// [[Rcpp::export]]
List boot_cb(NumericVector eval, NumericVector sam1, NumericVector sam0, 
                  int R, double eps) {
  std::sort(sam1.begin(), sam1.end());
  std::sort(sam0.begin(), sam0.end());
  int n1 = sam1.size();
  int n0 = sam0.size();
  int ne = eval.size();
  NumericVector F1 = edf(eval, sam1);
  NumericVector F0 = edf(eval, sam0);
  NumericVector s1(n1);
  NumericVector s0(n0);
  NumericMatrix boot(R, ne);
  NumericVector bF1(ne);
  NumericVector bF0(ne);
  NumericVector shift0(n0);
  NumericVector bshift0(n0);
  NumericVector dpro(ne);
  List ests = epsmax_cb(eval, sam1, sam0, eps);
  NumericMatrix emax = as<NumericMatrix>(ests["epsmaxL"]);
  //NumericMatrix emin = as<NumericMatrix>(ests["epsminU"]);
  NumericVector psiprimeL(ne);
  //NumericVector psiprimeU(ne);
  NumericVector simL_max(R);
  NumericVector simL_l2(R);
  //NumericVector simU(R);
  for (int i = 0; i < R; ++i) {
    s1 = sample(sam1, n1, true);
    s0 = sample(sam0, n0, true);
    std::sort(s1.begin(), s1.end());
    std::sort(s0.begin(), s0.end());
    bF1 = edf(eval, s1);
    for (int j = 0; j < ne; ++j) {
      for (int k = 0; k < n0; ++k) {
        shift0[k] = sam0[k] + eval[j];
        bshift0[k] = s0[k] + eval[j];
      }
      F0 = edf(eval, shift0);
      bF0 = edf(eval, bshift0);
      dpro = bF1 - F1 - bF0 + F0;
      psiprimeL[j] = max(as<NumericVector>(dpro[emax(_, j) == 1]));
      //psiprimeU[j] = min(as<NumericVector>(dpro[emin(_, j) == 1]));
    }
    simL_max[i] = max(abs(psiprimeL));
    simL_l2[i] = l2norm(eval, psiprimeL);
    //simU[i] = max(abs(psiprimeU));
  }
  return List::create(_["KS"] = simL_max, _["CvM"] = simL_l2);
}

