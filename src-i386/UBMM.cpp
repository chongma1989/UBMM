#ifndef __UBMM__
#define __UBMM__

// We can now use the BH package
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/special_functions/trigamma.hpp>

using namespace Rcpp;
using namespace boost::math;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::List UBMM(const NumericVector & x, const NumericVector & w, 
                const NumericVector & a, const double & precision,
                const int & iter=10000) {
  //Initialization
  NumericVector old_w=w;
  NumericVector old_a=a;
  NumericVector new_w(2),new_a(2);
  NumericVector logx=Rcpp::log(x),logxc=Rcpp::log(1-x);
  
  int itr=0;
  while(itr<iter){
    NumericVector old_dbeta=Rcpp::dbeta(x,old_a[0],old_a[1]);
    NumericVector weight=old_w[1]*old_dbeta/(old_w[0]+old_w[1]*old_dbeta);
    double sum_weight=Rcpp::sum(weight),mean_weight=Rcpp::mean(weight);
    
    //update weights
    new_w[0]=1-mean_weight;
    new_w[1]=1-new_w[0];
   
    //update parameters in Beta distribution
    new_a[0]=old_a[0]-((digamma(old_a[0]+old_a[1])-digamma(old_a[0]))*sum_weight+
      Rcpp::sum(weight*logx))/((trigamma(old_a[0]+old_a[1])-trigamma(old_a[0]))*sum_weight);
    
    new_a[1]=old_a[1]-((digamma(new_a[0]+old_a[1])-digamma(old_a[1]))*sum_weight+
      Rcpp::sum(weight*logxc))/((trigamma(new_a[0]+old_a[1])-trigamma(old_a[1]))*sum_weight);
    
    NumericVector new_dbeta=Rcpp::dbeta(x,new_a[0],new_a[1]);
    double loldsum=Rcpp::sum(Rcpp::log(old_w[0]+old_w[1]*old_dbeta));
    double lnewsum=Rcpp::sum(Rcpp::log(new_w[0]+old_w[1]*new_dbeta));
    
    if(std::fabs(loldsum-lnewsum)<precision){
      break;
    }else{
      old_a=new_a;
      old_w=new_w;
      itr=itr+1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Weight")=new_w,
                            Rcpp::Named("Beta Pars")=new_a,
                            Rcpp::Named("Iterations")=itr);
}


#endif 
