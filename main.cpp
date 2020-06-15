#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


mat matern_isotropic(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*Rf_gammafn(covparms(2)));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }

    
    
    // calculate covariances
    mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate distance
            double d = 0.0;
            for(int j=0; j<dim; j++){
                d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
            }
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance
                covmat(i2,i1) = normcon*
                    pow( d, covparms(2) )*Rf_bessel_k(d,covparms(2),1.0);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; }
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }
    }
    return covmat;
}



cube d_matern_isotropic(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    int n = locs.nrow();
    //double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*Rf_gammafn(covparms(2)));
    double eps = 1e-8;
    double normconeps =
        covparms(0)/(pow(2.0,covparms(2)+eps-1.0)*Rf_gammafn(covparms(2)+eps));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate derivatives
    cube dcovmat = cube(n,n, covparms.length(), fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        double cov;
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
            dcovmat(i1,i2,1) += 0.0;
            dcovmat(i1,i2,2) += 0.0;
        } else {
            cov = normcon*pow( d, covparms(2) )*Rf_bessel_k(d,covparms(2),1.0);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameter
            dcovmat(i1,i2,1) += normcon*pow(d,covparms(2))*
                Rf_bessel_k(d,covparms(2)-1.0,1.0)*d/covparms(1);
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,2) +=
                ( normconeps*pow(d,covparms(2)+eps)*Rf_bessel_k(d,covparms(2)+eps,1.0) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(3);
            dcovmat(i1,i2,3) += covparms(0);
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}



// [[Rcpp::export]]
List rlCompute(
     NumericVector covparms,
     NumericMatrix locs,
     IntegerMatrix NNarray,
     NumericVector y,
     NumericMatrix X,
     int r){
    
    long n = y.length();
    int m = NNarray.ncol();
    int p = X.ncol();
    int dim = locs.ncol();
    
    double rc  = 0.0;
    vec grad = zeros<vec>(covparms.length());
    mat info = zeros<mat>(covparms.length(), covparms.length());
    List ret;
    
    for(int i=0; i< n-r; i++){
       
       int bsize = std::min((r+i+1),m);

       NumericMatrix locsub(bsize, dim);
       vec ysub(bsize);
       mat X0( bsize, p );
        
       // form ysub,X0, and locs
     for(int j=bsize-1; j>=0; j--){
         ysub( (bsize-1) -j) = y( NNarray(r+i,j)-1 );
         for(int k=0;k<dim;k++){ locsub( (bsize-1) -j,k) = locs( NNarray(r+i,j)-1, k ); }
         for(int k=0;k<p;k++){ X0(bsize-1-j,k) = X( NNarray(r+i,j)-1, k ); }
     }
    
        
        
        // form covariance matrix and derivative
        mat covmat = matern_isotropic(covparms,locsub);
        cube dcovmat = d_matern_isotropic(covparms,locsub);
        
        // form peices for computing ith component of rl
        vec s = ysub.subvec(0,bsize-2);
        
        mat K =covmat.submat(0, 0, bsize-2, bsize-2);
        cube dK =dcovmat.subcube(0, 0,0,bsize-2, bsize-2,(covparms.length()-1));
        
        mat M = X0.submat(0,0,bsize-2,p-1);
        
        mat D =join_rows(join_cols(K, M.t()) ,join_cols(M, zeros<mat>(p,p)));
        
        vec k = (covmat.col(bsize-1)).subvec(0, bsize-2);
        vec b = join_cols(k,X0.row(bsize-1).t());
        
        vec lambda = solve(D,b).eval().head_rows(bsize-1);
        
        rowvec l = zeros<rowvec>(bsize);
        l.subvec(0, bsize-2)= -1*lambda.t();
        l.col(bsize-1) = 1;
        
        double w = ( l * ysub).eval()(0);
        double v = (l * covmat * l.t() ).eval()(0);
        
        // compute restricted likelihood
        rc +=   ((-0.5) *log(2 * M_PI)) + ((-0.5) *log(v)) + ((-0.5) *pow(w,2) *pow(v,-1)) ;
        
        vec gradTemp = zeros<vec>(covparms.length());
        mat infoTemp = zeros<mat>(covparms.length(), covparms.length());

        for(int k=0;k < (covparms.length()) ;k++){
            
            b = join_cols(dcovmat.slice(k).col(bsize-1).subvec(0, bsize-2) - dK.slice(k)*lambda , zeros<vec>(p));
            
            vec dlambda =(solve(D,b).eval()).head_rows(bsize-1);
            
            double dkwi = (  (-1*dlambda).t()  * s ).eval()(0);
                        
            double dkvi = ( l * dcovmat.slice(k) * l.t() ).eval()(0);
            
            gradTemp(k) =   (-0.5)* ((pow(v,-1) * dkvi)  + (2*w * pow(v,-1) * dkwi) - (pow(w,2) * pow(v,-2) * dkvi) );
                        
            for(int f=0;f< (covparms.length()) ;f++){
                double dfvi = ( l * dcovmat.slice(f) * l.t() ).eval()(0);
                infoTemp(k,f) = (0.5)*(pow(v,-2) * dkvi  * dfvi);
            }
        }
        
        info += infoTemp;
        grad += gradTemp;
    }
    ret = List::create( Named("loglik") = rc, Named("grad") = grad, Named("info") = info);
    return ret;
}
 
