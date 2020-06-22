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






void compute_pieces_2(
    NumericVector covparms,
    StringVector covfun_name,
    const NumericMatrix locs,
    IntegerMatrix NNarray,
    NumericVector& y,
    NumericMatrix X,
    mat* XSX,
    vec* ySX,
    double* ySy,
    double* wSw,
    double* logdet,
    cube* dXSX,
    mat* dySX,
    vec* dySy,
    vec* dwSw,
    vec* dlogdet,
    mat* ainfo,
    bool profbeta,
    bool grad_info,
    bool reml
){

    // data dimensions
    int n = y.length();
    int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    
    
    int r = 0;
    if(reml){
        r = p;
    }

    
    // loop over observations
    for(int i=0; i<n-r; i++){
    
        Rcpp::checkUserInterrupt();
        int bsize = std::min((r+i+1),m);

        // first, fill in ysub, locsub, and X0 in reverse order
        NumericMatrix locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0( bsize, p );
        for(int j=bsize-1; j>=0; j--){
            ysub( (bsize-1) -j) = y( NNarray(r+i,j)-1 );
            for(int k=0;k<dim;k++){ locsub( (bsize-1) -j,k) = locs( NNarray(r+i,j)-1, k ); }
            for(int k=0;k<p;k++){ X0(bsize-1-j,k) = X( NNarray(r+i,j)-1, k ); }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        
        arma::mat covmat = matern_isotropic( covparms, locsub );
        
        arma::cube dcovmat;
        if(grad_info){
            dcovmat = d_matern_isotropic( covparms, locsub );
        }
        
        arma::mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower" );
        
        // i1 is conditioning set, i2 is response
        //arma::span i1 = span(0,bsize-2);
        arma::span i2 = span(bsize-1,bsize-1);
        
        // get last row of cholmat
        arma::vec onevec = zeros(bsize);
        onevec(bsize-1) = 1.0;
        arma::vec choli2;
        
        if(grad_info){
            if(!reml){
            choli2 = solve( trimatu(cholmat.t()), onevec );
            }
        }
        
        bool cond = bsize > 1;
        //double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0;
        if(profbeta){
            LiX0 = solve( trimatl(cholmat), X0 );
        }
        arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        
        // loglik objects
        double v = 0.0;
        double w = 0.0;

        
        arma::vec s = arma::vec(bsize-1, fill::zeros);
        arma::mat K = mat( bsize-1, bsize-1,fill::zeros);
        arma::cube dK = cube(bsize-1, bsize-1,nparms, fill::zeros);
        arma::vec lambda = arma::vec(bsize-1, fill::zeros);
        arma::rowvec l = arma::rowvec(bsize, fill::zeros);
        arma::mat D = mat( (bsize-1) + p,(bsize-1) + p,fill::zeros );
        arma::mat M = mat(bsize-1,p,fill::zeros);
        
        if(reml){
            K = covmat.submat(0, 0, bsize-2, bsize-2);
            dK = dcovmat.subcube(0, 0,0,bsize-2, bsize-2,nparms-1);
            s = ysub.subvec(0,bsize-2);
            mat M = X0.submat(0,0,bsize-2,p-1);
            D = join_rows(join_cols(K, M.t()) ,join_cols(M, zeros<mat>(p,p)));
            vec k = (covmat.col(bsize-1)).subvec(0, bsize-2);
            lambda = solve(D,join_cols(k,X0.row(bsize-1).t())).eval().head_rows(bsize-1);
            l.subvec(0, bsize-2) = -1*lambda.t();
            l.col(bsize-1) = 1;
            w = ( l * ysub).eval()(0);
            v = (l * covmat * l.t() ).eval()(0);
            *logdet += std::log( v );
            *wSw +=   pow( w, 2 ) * pow(v,-1);
        }else{
            *logdet += 2.0*std::log( as_scalar(cholmat(i2,i2)) );
            *ySy +=    pow( as_scalar(Liy0(i2)), 2 );
        }
        
        
        if(profbeta){
            *XSX +=   LiX0.rows(i2).t() * LiX0.rows(i2);
            *ySX += ( Liy0(i2) * LiX0.rows(i2) ).t();
        }
        
        if( grad_info ){
        // gradient objects
        // LidSLi3 is last column of Li * (dS_j) * Lit for 1 parameter i
        // LidSLi2 stores these columns in a matrix for all parameters
        arma::mat LidSLi2(bsize,nparms);
        
        if(cond){ // if we condition on anything
            if(reml){
                 for(int j=0; j<nparms; j++){
                     vec b = join_cols(dcovmat.slice(j).col(bsize-1).subvec(0, bsize-2) - dK.slice(j)*lambda , zeros<vec>(p));
                     vec dlambda =(solve(D,b).eval()).head_rows(bsize-1);
                     double djw = (  (-1*dlambda).t()  * s ).eval()(0);;
                     double djv = ( l * dcovmat.slice(j) * l.t() ).eval()(0);;
                     (*dwSw)(j) += (2*w * pow(v,-1) * djw) - (pow(w,2) * pow(v,-2) * djv) ;
                     (*dlogdet)(j)+= (pow(v,-1) * djv);
                     
                     for(int f=0;f< nparms ;f++){
                         double dfv = ( l * dcovmat.slice(f) * l.t() ).eval()(0);
                         (*ainfo)(f,j) += (0.5)*(pow(v,-2) * djv  * dfv);
                     }
                }
            }
            
            else{
            for(int j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                arma::vec LidSLi3 = solve( trimatl(cholmat), dcovmat.slice(j) * choli2 );
                // store LiX0.t() * LidSLi3 and Liy0.t() * LidSLi3
                arma::vec v1 = LiX0.t() * LidSLi3;
                double s1 = as_scalar( Liy0.t() * LidSLi3 );
                // update all quantities
                // bottom-right corner gets double counted, so need to subtract it off
                (*dXSX).slice(j) += v1 * LiX0.rows(i2) + ( v1 * LiX0.rows(i2) ).t() -
                    as_scalar(LidSLi3(i2)) * ( LiX0.rows(i2).t() * LiX0.rows(i2) );
                (*dySy)(j) += as_scalar( 2.0 * s1 * Liy0(i2)  -
                    LidSLi3(i2) * Liy0(i2) * Liy0(i2) );
                (*dySX).col(j) += (  s1 * LiX0.rows(i2) + ( v1 * Liy0(i2) ).t() -
                    as_scalar( LidSLi3(i2) ) * LiX0.rows(i2) * as_scalar( Liy0(i2))).t();
                (*dlogdet)(j) += as_scalar( LidSLi3(i2) );
                // store last column of Li * (dS_j) * Lit
                LidSLi2.col(j) = LidSLi3;
            }

            // fisher information object
            // bottom right corner gets double counted, so subtract it off
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (*ainfo)(i,j) +=
                    1.0*accu( LidSLi2.col(i) % LidSLi2.col(j) ) -
                    0.5*accu( LidSLi2.rows(i2).col(j) %
                              LidSLi2.rows(i2).col(i) );
            }}
                }
        } else { // similar calculations, but for when there is no conditioning set
            for(int j=0; j<nparms; j++){
                arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                (*dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0;
                (*dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                (*dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                (*dlogdet)(j) += trace( LidSLi );
                LidSLi2.col(j) = LidSLi;
            }
            
            // fisher information object
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (*ainfo)(i,j) += 0.5*accu( LidSLi2.col(i) % LidSLi2.col(j) );
            }}

        }
        
        }

    }
}

void synthesize_2(
      NumericVector covparms,
      StringVector covfun_name,
      const NumericMatrix locs,
      IntegerMatrix NNarray,
      NumericVector& y,
      NumericMatrix X,
      NumericVector* ll,
      NumericVector* betahat,
      NumericVector* grad,
      NumericMatrix* info,
      NumericMatrix* betainfo,
      bool profbeta,
      bool grad_info,
      bool reml){

      // data dimensions
      int n = y.length();
      //int m = NNarray.ncol();
      int p = X.ncol();
      int nparms = covparms.length();
      //int dim = locs.ncol();
      
      
      
      // likelihood objects
      arma::mat XSX = arma::mat(p, p, fill::zeros);
      arma::vec ySX = arma::vec(p, fill::zeros);
      
      double ySy = 0.0;
      double wSw = 0.0;
      
      double logdet = 0.0;

      // gradient objects
      arma::cube dXSX = arma::cube(p,p,nparms,fill::zeros);
      arma::mat dySX = arma::mat(p, nparms, fill::zeros);
      
      arma::vec dwSw = arma::vec(nparms, fill::zeros);
      arma::vec dySy = arma::vec(nparms, fill::zeros);
      
      arma::vec dlogdet = arma::vec(nparms, fill::zeros);

      // fisher information
      arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);

      // this is where the big computation happens
      compute_pieces_2(
          covparms, covfun_name, locs, NNarray, y, X,
          &XSX, &ySX, &ySy, &wSw, &logdet, &dXSX, &dySX, &dySy, &dwSw, &dlogdet, &ainfo,
          profbeta, grad_info,reml
      );
          
      // synthesize everything and update loglik, grad, beta, betainfo, info
      
      // betahat and dbeta
      arma::vec abeta = arma::vec( p, fill::zeros );
      
      if(profbeta){ abeta = solve( XSX, ySX ); }
      for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); };

      arma::mat dbeta = arma::mat(p,nparms, fill::zeros);
      
      if( profbeta && grad_info){
      for(int j=0; j<nparms; j++){
          dbeta.col(j) = solve( XSX, dySX.col(j) - dXSX.slice(j) * abeta );
      }
      }
      
      if(reml){
          // loglikelihood
          (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + wSw );
      }else{
          double sig2 = ( ySy - 2.0*as_scalar( ySX.t() * abeta ) +
          as_scalar( abeta.t() * XSX * abeta ) )/n;
          // loglikelihood
          (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + n*sig2 );
      }
      
      
      if(profbeta){
      // betainfo
      for(int i=0; i<p; i++){ for(int j=0; j<i+1; j++){
          (*betainfo)(i,j) = XSX(i,j);
          (*betainfo)(j,i) = XSX(j,i);
      }}
      }

      if(grad_info){
          if(reml){
              //gradient
              for(int j=0; j<nparms; j++){
                  (*grad)(j) = 0.0;
                  (*grad)(j) -= 0.5*dlogdet(j);
                  (*grad)(j) -= 0.5*dwSw(j);
                  // fisher information
                  for(int i=0; i<nparms; i++){
                      for(int j=0; j<i+1; j++){
                      (*info)(i,j) = ainfo(i,j);
                      (*info)(j,i) = (*info)(i,j);
                      }
                  }
              }
          }
          else{
              for(int j=0; j<nparms; j++){
                  (*grad)(j) = 0.0;
                  (*grad)(j) -= 0.5*dlogdet(j);
                  (*grad)(j) += 0.5*dySy(j);
                  (*grad)(j) -= 1.0*as_scalar( abeta.t() * dySX.col(j) );
                  (*grad)(j) += 1.0*as_scalar( ySX.t() * dbeta.col(j) );
                  (*grad)(j) += 0.5*as_scalar( abeta.t() * dXSX.slice(j) * abeta );
                  (*grad)(j) -= 1.0*as_scalar( abeta.t() * XSX * dbeta.col(j) );
              }
              // fisher information
              for(int i=0; i<nparms; i++){
                  for(int j=0; j<nparms; j++){
                 (*info)(i,j) = ainfo(i,j);
                 (*info)(j,i) = (*info)(i,j);
                  }
              }
          }
      }
  }


// [[Rcpp::export]]
List vecchia_reml_loglik_grad_info(
    NumericVector covparms,
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    synthesize_2(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo, false, true, true
    );
    
    List ret = List::create( Named("loglik") = ll,
        Named("grad") = grad, Named("info") = info );
    return ret;
        
}


// [[Rcpp::export]]
List vecchia_profbeta_loglik_grad_info_2(
    NumericVector covparms,
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    synthesize_2(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo, true, true,false
    );
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("grad") = grad, Named("info") = info, Named("betainfo") = betainfo );
    return ret;
        
}





