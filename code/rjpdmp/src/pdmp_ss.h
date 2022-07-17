#ifndef PDMP_SS_H
#define PDMP_SS_H

#include "inversion.h"
#include "modelprior.h"
#include "likelyhood.h"

arma::vec zigzag_MP_GaussIID_Mu(arma::vec theta, const arma::vec& mu,double rjval, double ppi, double prior_sigma2) {
  arma::uvec inds_off = arma::find(arma::abs(theta) > 1e-10);
  arma::uvec inds_on = arma::find(arma::abs(theta) < 1e-10);
  int p = inds_on.size();
  arma::vec times(p);
  double alpha;

  for(int i =0; i < p; i++){
    alpha = rjval*ppi/(1-ppi)/sqrt(prior_sigma2*2.0*arma::datum::pi);
    alpha = alpha*exp(-0.5*std::pow(mu(i),2.0));
    times(i) = R::rexp(1)/alpha;
  }
  return(times);
}

arma::vec bps_MP_GaussIID_N_Mu(arma::vec theta,const arma::vec& mu, double rjval, double ppi, double prior_sigma2) {
  arma::uvec inds_off = arma::find(arma::abs(theta) > 1e-10);
  arma::uvec inds_on = arma::find(arma::abs(theta) < 1e-10);

  int num_off = inds_off.size();
  int num_on = theta.size() - num_off;
  double Jac = 1.0;

  if( num_off > 1e-10 ){
    Jac = 2.0/sqrt(2.0*arma::datum::pi);
  }

  arma::vec times(num_on);
  double alpha;
  for(int i=0; i < num_on; i++){
    alpha = -0.5*log(2*arma::datum::pi*prior_sigma2) -0.5*std::pow(mu(i),2.0);
    alpha += log(rjval) + log(ppi)-log(1-ppi) + log(Jac);
    times(i) = R::rexp(1)/exp(alpha);
  }
  return(times);
}
double logarea_s(double d){
  return(2 + (0.5*d+0.5)*log(arma::datum::pi) - lgamma(0.5*d+0.5));
}

arma::vec bps_MP_GaussIID_S_Mu(arma::vec theta,const arma::vec& mu, double rjval, double ppi, double prior_sigma2) {
  arma::uvec inds_off = arma::find(arma::abs(theta) > 1e-10);
  arma::uvec inds_on = arma::find(arma::abs(theta) < 1e-10);

  int num_off = inds_off.size();
  int num_on = theta.size() - num_off;
  double lJac = 0.0;

  if( num_off > 1e-10 ){
    lJac = logarea_s(theta.size() - num_on -1) - logarea_s(theta.size() - num_on);
    lJac = lJac + log(2.0) -log(theta.size() - num_on);

  }
  arma::vec times(num_on);

  double alpha;
  for(int i=0; i < num_on; i++){
    alpha = -0.5*log(2*arma::datum::pi*prior_sigma2) -0.5*std::pow(mu(i),2.0);
    alpha += log(rjval) + log(ppi)-log(1-ppi) + lJac;
    times(i) = R::rexp(1)/exp(alpha);
  }
  return(times);
}


//' zigzag_ss
//'
//' Applies the Reversible Jump ZigZag Sampler to a Dirac spike and slab distribution.
//'
//' @param maxTime Maximum runtime (in Seconds) of the algorithm; will terminate the code after a given computation time or nmax iterations of the algorithm is reached.
//' @param sigma2 Double for the prior variance for included variables.
//' @param mu Double for the prior variance for included variables.
//' @param x0 Initial position of the regression parameter
//' @param theta0 Initial velocity for the sampler (Default has 1s on all components). This should be chosen with unit velocities on each component (regardless of sign).
//' @param rj_val Reversible jump parameter for the PDMP method. This value is fixed over all models and is interpreted as the probability to jump to a reduced model when a parameter hits zero.
//' @param ppi Double for the prior probability of inclusion (ppi) for each parameter.
//' @param nmax Maximum number of iterations (simulated events) of the algorithm; will stop the algorithm when this number of iterations of the method have occured. Default value is 1e6, lower values should be chosen for memory constraints if less iterations are desired.
//' @param burn Optional number of iterations to use for burnin. These are not stored so can be useful in memory intensive problems.
//' @return Returns a list with the following objects:
//' @return \code{times}: Vector of event times where ZigZag process switchs velocity or jumps models.
//' @return \code{positions}: Matrix of positions at which event times occur, these are not samples from the PDMP.
//' @return \code{theta}: Matrix of new velocities at event times.
//' @examples
//'
//' p <- 100
//' sig <- 1
//' ppi <- 0.5
//'
//' @export
// [[Rcpp::export]]
List zigzag_ss(double maxTime, const arma::vec& mu,
                    double sigma2, arma::vec x0, arma::vec theta0,
                    double rj_val = 0.6, double ppi=0.5, int nmax = 1e6, int burn = -1){
  int mini = 1, p = x0.size(), nEvent= 1;
  double eps = 1e-10, t = 0, upper, val, tau_val;

  arma::mat ab_vals(p,2), sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p), x = x0, grad_val(1), clock_times(nmax);

  taus.zeros(); ab_vals.zeros(); sk_times.zeros(), clock_times.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }

  arma::vec grad_vals(p), a_vals(p), b_vals(p);
  b_vals.zeros(); a_vals.zeros();

  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);
  arma::uvec miniu(1);

  grad_vals.elem(inds_off_hp) = (x.elem(inds_off_hp) - mu(inds_off_hp))/sigma2;
  a_vals(inds_off_hp) = theta(inds_off_hp)%grad_vals.elem(inds_off_hp);
  b_vals(inds_off_hp) = theta(inds_off_hp)%theta(inds_off_hp)/sigma2;

  taus(inds_off_hp) = get_time_off_hp(x(inds_off_hp),theta(inds_off_hp),
       a_vals(inds_off_hp), b_vals(inds_off_hp));

  taus(inds_on_hp) = zigzag_MP_GaussIID_Mu(theta,mu, rj_val, ppi, sigma2);

  while( nEvent < nmax + burn ){
    mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    a_vals += b_vals*tau_val;

    if( std::abs(x(mini)) < eps){
      x(mini) = 0;
      if( std::abs(theta(mini)) < eps ){
        // Exiting the HP
        theta[mini] = 2.0*(R::runif(0,1) < 0.5) - 1;
        if( nEvent >= burn){
          clock_times(nEvent-burn) = timer.toc();
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      } else {
        // Enter the HP
        // w.p. rj enter the HP
        if(R::runif(0,1) < rj_val){
          theta[mini] = 0;
          if( nEvent >= burn){
            clock_times(nEvent-burn) = timer.toc();
            sk_times(nEvent-burn) = t;
            sk_points.col(nEvent-burn) = x;
            sk_theta.col(nEvent-burn) = theta;
          }
          nEvent++;
        }
      }
      // inds_off_hp = arma::find(arma::abs(theta) > eps);
      // inds_on_hp = arma::find(arma::abs(theta) < eps);
      //
      // grad_vals.elem(inds_off_hp) = (x.elem(inds_off_hp) - mu(inds_off_hp))/sigma2;
      // a_vals(inds_off_hp) = theta(inds_off_hp)%grad_vals.elem(inds_off_hp);
      // b_vals(inds_off_hp) = theta(inds_off_hp)%theta(inds_off_hp)/sigma2;

    } else {

      upper = std::max(0.0, a_vals(mini));
      miniu(0) = mini;
      grad_val = (x(mini) - mu(mini))/sigma2;

      val = theta(mini)*grad_val(0);
      if( R::runif(0,1) < val/upper){
        theta(mini) = -theta(mini);
        if( nEvent >= burn){
          clock_times(nEvent-burn) = timer.toc();
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
        b_vals(mini) = theta(mini)*theta(mini)/sigma2;
      }
      a_vals(mini) = val;
    }
    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime" << std::endl;
      } else {
        clock_times(nEvent-burn) = timer.toc();
        sk_points.shed_cols(nEvent-burn, nmax-1);
        sk_theta.shed_cols(nEvent-burn, nmax-1);
        sk_times.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }

    inds_off_hp = arma::find(arma::abs(theta) > eps);
    inds_on_hp = arma::find(arma::abs(theta) < eps);

    grad_vals.elem(inds_off_hp) = (x.elem(inds_off_hp) - mu(inds_off_hp))/sigma2;
    a_vals(inds_off_hp) = theta(inds_off_hp)%grad_vals.elem(inds_off_hp);
    b_vals(inds_off_hp) = theta(inds_off_hp)%theta(inds_off_hp)/sigma2;

    taus(inds_off_hp) = t + get_time_off_hp(x(inds_off_hp),theta(inds_off_hp),
         a_vals(inds_off_hp), b_vals(inds_off_hp));

    taus(inds_on_hp) = t + zigzag_MP_GaussIID_Mu(theta,mu, rj_val, ppi, sigma2);

  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  ret["clock"] = clock_times ;
  return(ret) ;
}

// [[Rcpp::export]]
List bps_n_ss(double maxTime, const arma::vec& mu,
                 double sigma2, arma::vec x0, arma::vec theta0,
                 double ref0 = 0.1, double rj_val = 0.6,
                 double ppi=0.5, int nmax = 1e6, int burn = -1, bool ref_auto = true){

  int p = x0.size(), nEvent= 1;
  double eps = 1e-10, epsbig = 0.5, t = 0, upper, val, tau_val, alpha;

  arma::mat sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p+2), x = x0, en_variables(p), ab_vals(2);

  taus.zeros(); ab_vals.zeros(); sk_times.zeros(), en_variables.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }
  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  // Assume linear bounds
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);

  en_variables.elem(inds_off_hp).ones();
  arma::vec grad_vals(p);
  grad_vals.zeros();

  grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;
  ab_vals(0) = arma::dot(grad_vals(inds_off_hp), theta(inds_off_hp));
  ab_vals(1) = arma::dot(theta(inds_off_hp), theta(inds_off_hp))/sigma2;

  taus(inds_off_hp) = get_hit_times(x(inds_off_hp),theta(inds_off_hp));
  taus(inds_on_hp) = bps_MP_GaussIID_N_Mu(theta, mu, rj_val, ppi, sigma2);
  taus(p) = linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));

  double ref = ref0;
  if((ref_auto) & (inds_off_hp.size() > 0.0)){
    ref = ref0*sqrt(inds_off_hp.size());
  }
  taus(p+1) = R::rexp(1)/ref;

  // Rcout << " as ";
  while( nEvent < nmax + burn){
    int mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    ab_vals(0) += ab_vals(1)*tau_val;

    if( mini == p + 1){
      // Refresh
      theta.elem(arma::find(en_variables > epsbig)).randn();
      if( nEvent >= burn){
        sk_times(nEvent-burn) = t;
        sk_points.col(nEvent-burn) = x;
        sk_theta.col(nEvent-burn) = theta;
      }
      nEvent++;
      grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;
    }
    if( mini == p){
      // Flip operator
      upper = std::max(0.0, ab_vals(0));
      grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;

      val = arma::dot(grad_vals(inds_off_hp),theta(inds_off_hp));
      if( R::runif(0,1) < val/upper){

        theta.elem(inds_off_hp) -=
          2*arma::dot(theta.elem(inds_off_hp), grad_vals.elem(inds_off_hp))/
            arma::dot(grad_vals.elem(inds_off_hp),grad_vals.elem(inds_off_hp))*grad_vals.elem(inds_off_hp);

        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }
    }
    if( mini < p) {
      // Rcout << "3";
      x(mini) = 0;
      inds_off_hp = arma::find(en_variables > epsbig);
      if( en_variables(mini) > epsbig ){
        // w.p. rj enter HP
        if(R::runif(0,1) < rj_val){
          theta(mini) = 0;
          en_variables(mini) = 0;
          if( nEvent >= burn){
            sk_times(nEvent-burn) = t;
            sk_points.col(nEvent-burn) = x;
            sk_theta.col(nEvent-burn) = theta;
          }
          nEvent++;
        }
      } else {
        // if mini is on HP then exit the HP
        if( inds_off_hp.size()  < epsbig ){
          alpha = 2.0*(R::runif(0,1) < 0.5) - 1;
        } else {
          alpha = MH_alpha_N(inds_off_hp.size());
        }
        theta(mini) = alpha;

        en_variables(mini) = 1;
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }

      inds_off_hp = arma::find(en_variables > epsbig);
      inds_on_hp = arma::find(en_variables < epsbig);

      grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;
    }

    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime" << std::endl;
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        sk_theta.shed_cols(nEvent-burn, nmax-1);
        sk_times.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
    if(arma::sum(en_variables) > 0.0){
      ab_vals(0) = arma::dot(grad_vals(inds_off_hp), theta(inds_off_hp));
      ab_vals(1) = arma::dot(theta(inds_off_hp), theta(inds_off_hp))/sigma2;
      taus.elem(inds_off_hp) = t + get_hit_times(x.elem(inds_off_hp),theta.elem(inds_off_hp));
      taus(p) = t + linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));
    } else {
      taus(p) = t + 1e10;
    }
    if(arma::sum(en_variables) < p){
      taus.elem(inds_on_hp) = t + bps_MP_GaussIID_N_Mu(theta, mu, rj_val, ppi, sigma2);
    }
    if((ref_auto) & (arma::sum(en_variables) > 0.0)){
      ref = ref0*sqrt(arma::sum(en_variables));
    }
    taus(p+1) = t + R::rexp(1)/ref;


  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  return(ret) ;
}

// [[Rcpp::export]]
List bps_s_ss(double maxTime, const arma::vec& mu,
              double sigma2, arma::vec x0, arma::vec theta0,
              double ref0 = 0.1, double rj_val = 0.6,
              double ppi=0.5, int nmax = 1e6, int burn = -1, bool ref_auto = true){

  int p = x0.size(), nEvent= 1;
  double eps = 1e-10, epsbig = 0.5, t = 0, upper, val, tau_val, alpha;

  arma::mat sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p+2), x = x0, en_variables(p), ab_vals(2);

  taus.zeros(); ab_vals.zeros(); sk_times.zeros(), en_variables.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }
  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  // Assume linear bounds
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);

  en_variables.elem(inds_off_hp).ones();
  if(sum(en_variables) > epsbig ){
    theta = theta/sqrt(arma::dot(theta,theta));
  }
  arma::vec grad_vals(p);
  grad_vals.zeros();

  grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;
  ab_vals(0) = arma::dot(grad_vals(inds_off_hp), theta(inds_off_hp));
  ab_vals(1) = arma::dot(theta(inds_off_hp), theta(inds_off_hp))/sigma2;

  taus(inds_off_hp) = get_hit_times(x(inds_off_hp),theta(inds_off_hp));
  taus(inds_on_hp) = bps_MP_GaussIID_S_Mu(theta, mu, rj_val, ppi, sigma2);
  taus(p) = linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));

  double ref = ref0;
  if((ref_auto) & (inds_off_hp.size() > 0.0)){
    ref = ref0*sqrt(inds_off_hp.size());
  }
  taus(p+1) = R::rexp(1)/ref;

  // Rcout << " as ";
  while( nEvent < nmax + burn){
    int mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    ab_vals(0) += ab_vals(1)*tau_val;

    if( mini == p + 1){
      // Refresh
      inds_off_hp = arma::find(en_variables > epsbig);
      theta.elem(inds_off_hp).randn();
      if(sum(en_variables) > epsbig ){
        theta(inds_off_hp) = theta(inds_off_hp)/sqrt(arma::dot(theta(inds_off_hp),theta(inds_off_hp)));
      }
      if( nEvent >= burn){
        sk_times(nEvent-burn) = t;
        sk_points.col(nEvent-burn) = x;
        sk_theta.col(nEvent-burn) = theta;
      }
      nEvent++;
      grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;
    }
    if( mini == p){
      // Flip operator
      upper = std::max(0.0, ab_vals(0));
      grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;

      val = arma::dot(grad_vals(inds_off_hp),theta(inds_off_hp));
      if( R::runif(0,1) < val/upper){

        theta.elem(inds_off_hp) -=
          2*arma::dot(theta.elem(inds_off_hp), grad_vals.elem(inds_off_hp))/
            arma::dot(grad_vals.elem(inds_off_hp),grad_vals.elem(inds_off_hp))*grad_vals.elem(inds_off_hp);

        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }
    }
    if( mini < p) {
      // Rcout << "3";
      x(mini) = 0;
      inds_off_hp = arma::find(en_variables > epsbig);
      if( en_variables(mini) > epsbig ){
        // w.p. rj enter HP
        if(R::runif(0,1) < rj_val){
          theta(mini) = 0;
          en_variables(mini) = 0;
          inds_off_hp = arma::find(en_variables > epsbig);
          if(sum(en_variables) > epsbig ){
            theta(inds_off_hp) = theta(inds_off_hp)/sqrt(arma::dot(theta(inds_off_hp),theta(inds_off_hp)));
          }
          if( nEvent >= burn){
            sk_times(nEvent-burn) = t;
            sk_points.col(nEvent-burn) = x;
            sk_theta.col(nEvent-burn) = theta;
          }
          nEvent++;
        }
      } else {
        // if mini is on HP then exit the HP
        if( inds_off_hp.size()  < epsbig ){
          alpha = 2.0*(R::runif(0,1) < 0.5) - 1;
        } else {
          alpha = sim_alpha(inds_off_hp.size());
        }
        theta = sqrt(1-std::pow(alpha,2))*theta;
        theta[mini] = alpha;
        en_variables(mini) = 1;

        inds_off_hp = arma::find(en_variables > epsbig);
        theta(inds_off_hp) = theta(inds_off_hp)/sqrt(arma::dot(theta(inds_off_hp),theta(inds_off_hp)));


        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }

      inds_off_hp = arma::find(en_variables > epsbig);
      inds_on_hp = arma::find(en_variables < epsbig);

      grad_vals(inds_off_hp) = (x(inds_off_hp) - mu(inds_off_hp))/sigma2;
    }

    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime" << std::endl;
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        sk_theta.shed_cols(nEvent-burn, nmax-1);
        sk_times.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
    if(arma::sum(en_variables) > 0.0){
      ab_vals(0) = arma::dot(grad_vals(inds_off_hp), theta(inds_off_hp));
      ab_vals(1) = arma::dot(theta(inds_off_hp), theta(inds_off_hp))/sigma2;
      taus.elem(inds_off_hp) = t + get_hit_times(x.elem(inds_off_hp),theta.elem(inds_off_hp));
      taus(p) = t + linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));
    } else {
      taus(p) = t + 1e10;
    }
    if(arma::sum(en_variables) < p){
      taus.elem(inds_on_hp) = t + bps_MP_GaussIID_S_Mu(theta, mu, rj_val, ppi, sigma2);
    }
    if((ref_auto) & (arma::sum(en_variables) > 0.0)){
      ref = ref0*sqrt(arma::sum(en_variables));
    }
    taus(p+1) = t + R::rexp(1)/ref;


  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  return(ret) ;
}

#endif
