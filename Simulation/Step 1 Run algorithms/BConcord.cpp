#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
List BCONCORDSS(mat S, int nmc, int burnin, int n, mat omega_start){

  int p = S.n_cols;
  mat Omega = omega_start;
  
  double pi = 3.141592;
  cube Bs = zeros(p,p,2); // Bs is an array of the votes given to zero and non-zero edges
  vec P = zeros(2);
  vec myvec = zeros(2);
  for(int t=0; t<2;t++)
    myvec(t) = t;
  double a = 0;
  double b = 0;
  double bd = 0;
  double lambda = 0;
  int rand = 1;
  
  double r = 1e-4;
  double s = 1e-8;

  for(int it = 0; it < (nmc + burnin); it++){
    
    for(int j=0; j<p-1; j++){
      for(int k=j+1; k<p; k++){
        
        a = S(j,j) + S(k,k);
        b = - sum(Omega.col(j)%S.col(k) + Omega.col(k)%S.col(j)) + Omega(j,k)*a;
        
        lambda = R::rgamma(r + 0.5, 1/(0.5*Omega(j,k)*Omega(j,k) + s));

        P(1) = sqrt(2*pi/(n*a+lambda)) * exp((n*b)*(n*b)/(2*(n*a+lambda))); 
        P(0) = 1;
        
        if(P.has_inf()==true){
          int rand = 1;
          Omega(j,k) = Omega(k,j) = R::rnorm(n*b/(n*a+lambda), 1/(n*a + lambda));
          
          if((it > burnin))
            Bs(j,k,rand) += 1;
            
        }  else {
          
          
          P = P/sum(P);
          
          rand = (Rcpp::RcppArmadillo::sample(myvec,1,false, P))(0);

          if(rand == 1){
            Omega(j,k) = Omega(k,j) = R::rnorm(n*b/(n*a+lambda), 1/(n*a + lambda));
          } else {
            Omega(j,k) = Omega(k,j) = 0;
          }
     
          if((it > burnin))
            Bs(j,k,rand) += 1; 
            
        }

        
      }
      
      bd = sum(Omega.col(j)%S.col(j)) - Omega(j,j)*S(j,j) ;
      Omega(j,j) = (sqrt(bd*bd + 4*S(j,j)) - bd)/(2*S(j,j));
    }
    bd = sum(Omega.col(p-1)%S.col(p-1)) - Omega(p-1,p-1)*S(p-1,p-1);
    Omega(p-1,p-1) = (sqrt(bd*bd + 4*S(p-1,p-1)) - bd)/(2*S(p-1,p-1));
    
    
    
  }
  
  mat edge_incl_prob = zeros(p,p);
    
  for(int j=0; j<p-1; j++){
    for(int k=j+1; k<p; k++){
      edge_incl_prob(j,k) = Bs(j,k,1)/nmc;
      edge_incl_prob(k,j) = edge_incl_prob(j,k);
    }
  }
  
  return(List::create(Named("p_links") = edge_incl_prob, Named("last_omega") = Omega));
}