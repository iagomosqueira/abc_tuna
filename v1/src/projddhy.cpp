//////////////////////////////////////////////////////
// projdd.cpp ////////////////////////////////////////
//////////////////////////////////////////////////////
// R. Hillary CSIRO 2021 /////////////////////////////
//////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <R.h>
#include "assert.h"

using namespace Rcpp;

/*
 * Project forward in time for fixed exploitation rate
 * inputs:
 * 1.  Ninit: initial population
 * 2.  B0
 * 3.  steepness
 * 4.  sigmar
 * 5.  rho
 * 6.  M
 * 7.  weight-at-age
 * 8.  maturity-at-age
 * 9.  selectivity-at-age
 * 10. fleet F share
 * 11. fixed explotation rate
 * 12. dimensions
 */
//[[Rcpp::export]]
RcppExport SEXP projddhy(SEXP Ninit_,SEXP B0_,SEXP hh_,SEXP SPR_,SEXP sigmar_,SEXP rho_,SEXP M_,SEXP wt_,SEXP mat_,SEXP sel_,SEXP fshare_,SEXP xiy_,SEXP dm_)
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  int ny = dm[0];
  int na = dm[1];
  int nf = dm[2];
  int a,f,i,y;

  // read ins

  NumericVector Ninit = as<Rcpp::NumericVector>(Ninit_);
  NumericVector M = as<Rcpp::NumericVector>(M_);
  NumericVector wt = as<Rcpp::NumericVector>(wt_); 
  NumericVector mat = as<Rcpp::NumericVector>(mat_);
  NumericMatrix sel = as<Rcpp::NumericMatrix>(sel_); 
  NumericVector fshare = as<Rcpp::NumericVector>(fshare_); 
  NumericVector xiy = as<Rcpp::NumericVector>(xiy_);
  double B0 = as<double>(B0_);
  double hh = as<double>(hh_);
  double SPR = as<double>(SPR_);
  double sigmar = as<double>(sigmar_);
  double rho = as<double>(rho_); 
  double alp = 4.*hh/(SPR*(1.-hh));
  double bet = (5.*hh-1.)/(B0*(1.-hh));
  double eps1,eps2,xisum;

  // memory allocation

  NumericMatrix N(na,ny);
  NumericMatrix N0(na,ny);
  NumericVector SSB(ny);
  NumericVector SSB0(ny); 

  // set up the fixed exploitation rates

  //for(f=0;f<nf;f++)
    //for(a=0;a<na;a++) xi(a,f) = xifix*fshare(f)*sel(f,a);

  // init popn

  for(a=0;a<na;a++) {
    
    N(a,0) = Ninit(a);
    N0(a,0) = Ninit(a);

  }

  SSB(0) = SSB0(0) = 0.;
  for(a=0;a<na;a++) {
    
    SSB(0) += N(a,0)*wt(a)*mat(a);
    SSB0(0) += N0(a,0)*wt(a)*mat(a); 

  }

  // loop through the years

  eps1 = R::rnorm(-sigmar*sigmar/2.,sigmar);
  for(y=1;y<ny;y++) {

    // recruitment

    if(y == 1) {
      
      N(0,y) = alp*SSB(y-1)/(1.+bet*SSB(y-1))*exp(eps1);
      N0(0,y) = N(0,y);

    } else {

      eps2 = rho*eps1+R::rnorm(-sigmar*sigmar/2.,sqrt(1.-rho*rho)*sigmar);
      N(0,y) = alp*SSB(y-1)/(1.+bet*SSB(y-1))*exp(eps2); 
      N0(0,y) = N(0,y); 

    }

    eps1 = eps2;

    // rest of the ages

    for(a=1;a<na;a++) {

      for(xisum=0.,f=0;f<nf;f++) xisum += xiy(y-1)*fshare(f)*sel(a-1,f);
      N(a,y) = N(a-1,y-1)*exp(-M(a-1))*(1.-xisum);
      N0(a,y) = N0(a-1,y-1)*exp(-M(a-1)); 

    }

    // plus group

    for(xisum=0.,f=0;f<nf;f++) xisum += xiy(y-1)*fshare(f)*sel(na-1,f); 
    N(na-1,y) += N(na-1,y-1)*exp(-M(na-1))*(1.-xisum);
    N0(na-1,y) += N0(na-1,y-1)*exp(-M(na-1));
    
    // SSB

    SSB(y) = SSB0(y) = 0.;
    for(a=0;a<na;a++) {
      
      SSB(y) += N(a,y)*wt(a)*mat(a); 
      SSB0(y) += N0(a,y)*wt(a)*mat(a); 

    }
  }

  List res = Rcpp::List::create(Named("N")=N,Named("N0")=N0);

  return Rcpp::wrap(res);

}
