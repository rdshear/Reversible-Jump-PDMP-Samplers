#ifndef HMC_LOGIT_H
#define HMC_LOGIT_H

// Negative log posterior
// [[Rcpp::export]]
double ll_logit(const arma::mat& dataX, const arma::vec& datay, arma::vec x, arma::vec gamma, double prior_sigma2){
  double s = sum(gamma);
  if( s < 1 ){
    return(sum( log(1+exp(datay*0)) - datay*0));
    //return(0.0);
  }
  arma::uvec nzero = find(arma::abs(gamma) > 1e-10);
  arma::mat dataX_sub = dataX.cols(nzero);
  // double log_post = sum( datay%log(1+exp(-dataX_sub*x(nzero))) +(1-datay)%(dataX_sub*x(nzero) + log(1+ exp(-dataX_sub*x(nzero))))) +
  //   0.5*arma::dot(x(nzero),x(nzero))/prior_sigma2 + (s/2)*log(2*arma::datum::pi);

  double log_post = sum( log(1+exp(dataX_sub*x(nzero))) - datay%(dataX_sub*x(nzero))) +
    0.5*arma::dot(x(nzero),x(nzero))/prior_sigma2 + (s/2)*log(2*arma::datum::pi*prior_sigma2);
  return(log_post);
}
double q_normal(arma::vec xv, double rw_sig2){
  if(xv.size() < 1){
    return(1);
  } else {
    return(exp(-0.5*arma::dot(xv,xv/rw_sig2) - 0.5*xv.size()*log(2*arma::datum::pi*rw_sig2) ));
  }
}

// [[Rcpp::export]]
double calc_gamma(arma::vec gamma, arma::vec gamma_change){

  arma::uvec gamma_on = find(arma::abs(gamma) > 1e-10);
  arma::uvec gamma_change_on = find(arma::abs(gamma_change) > 1e-10);
  int p = gamma.size();

  if(gamma_on.size() < 1){
    return(0.5);
  } else {
    if(gamma_on.size() == p){
      return(0.5);
    } else {
      return(1.0/3.0);
    }
  }
}

// [[Rcpp::export]]
List hmc_logit_rj(double maxTime, const arma::mat& X, const arma::vec& y, arma::vec x0,
               double epsilon = 0.1, double stoch_time = 1, int nmax = 10^6,
               int burn = -1, int thin = 1, double prior_sigma2 = 1, double ppi = 0.5, double rw_sig2 = 0.5){

  arma::wall_clock timer;
  timer.tic();

  int p = x0.size(), nEvent = 1, rj_birth = 1/3, rj_death = 1/3;
  int L = std::ceil(stoch_time/epsilon), g_i = 1, move = 1;
  double eps = 1e-10, acc_prob, Ux, Ux_prop, Uz, Uz_prop, prior_diff_gamma;
  double s_gamma, ratio_val, p_gamma_change, q_norm_rj, adj_alpha;

  arma::mat sk_points(p,nmax);
  arma::vec acc_count(nmax), x = x0, x_prop(p), z_vals(p), z_dash(p), gamma(p), gamma_prop(p);
  gamma.zeros();

  if( burn < 0){
    burn = 0;   sk_points.col(0) = x0;
  }
  arma::uvec inds_off_hp = arma::find(arma::abs(x) > eps);
  gamma(inds_off_hp).ones();
  Ux = ll_logit(X, y, x, gamma, prior_sigma2);
  arma::vec x_grad(p);

  arma::uvec inds_on_hp = arma::find(arma::abs(gamma) < eps);

  while( nEvent < nmax + burn){

    for( int it = 0; it < thin; it++){
      // Determine move: 1 = within, 2 = birth, 3 = death
      if(inds_off_hp.size() == p){
        move = 3;
        if(R::runif(0,1) < 0.5){
          move = 1;
        }
      } else {
        if(inds_off_hp.size() < eps){
          move = 2;
          if(R::runif(0,1) < 0.5){
            move = 1;
          }
        } else {
          move = floor(3.0*R::runif(0,1) + 1.0);
        }
      }

      // Within Model
      if(move == 1){
        if(sum(gamma) > eps){
          z_vals = rnorm(sum(gamma));
          x_prop = x;
          Uz = arma::dot(z_vals,z_vals)/2.0;
          z_dash = z_vals;
          x_grad.zeros();
          for( int l_val = 0; l_val < L; l_val++){
            x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
            z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
            x_prop(inds_off_hp) += epsilon * z_dash;

            x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
            z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
          }
          z_dash = -z_dash;
          Uz_prop = arma::dot(z_dash,z_dash)/2.0;
          Ux_prop = ll_logit(X, y, x_prop, gamma, prior_sigma2);

          acc_prob = std::min(1.0, exp(Ux - Ux_prop + Uz - Uz_prop));
          if( R::runif(0,1) < acc_prob ){
            x = x_prop;
            Ux = Ux_prop;
          }
        }
      } else {
        x_prop = x;
        gamma_prop = gamma;

        // Birth
        if(move == 2){
          g_i = floor(inds_on_hp.size()*R::runif(0,1));
          gamma_prop[inds_on_hp[g_i]] = 1.0;
          x_prop[inds_on_hp[g_i]] = R::rnorm(0, sqrt(rw_sig2));
          prior_diff_gamma = log(1-ppi) - log(ppi);

          q_norm_rj = -R::dnorm(x_prop[inds_on_hp[g_i]], 0.0, sqrt(rw_sig2), true);
        }
        // Death
        if(move == 3){
          g_i = floor(inds_off_hp.size()*R::runif(0,1));
          gamma_prop[inds_off_hp[g_i]] = 0.0;
          prior_diff_gamma = log(ppi) - log(1-ppi);
          q_norm_rj = R::dnorm(x_prop[inds_off_hp[g_i]], 0.0, sqrt(rw_sig2), true);
        }
        ratio_val = -ll_logit(X, y, x_prop, gamma_prop, prior_sigma2)+ll_logit(X, y, x, gamma, prior_sigma2);
        adj_alpha = calc_gamma(gamma, gamma_prop)/calc_gamma(gamma_prop, gamma);
        p_gamma_change = exp(ratio_val-prior_diff_gamma + q_norm_rj + log(adj_alpha));

        if(R::runif(0,1) < p_gamma_change){
          x = x_prop;
          gamma = gamma_prop;
          Ux = ll_logit(X, y, x, gamma, prior_sigma2);//ratio_val + Ux;

          inds_off_hp = arma::find(arma::abs(gamma) >= eps);
          inds_on_hp = arma::find(arma::abs(gamma) < eps);
        }
      }
    }

    if( nEvent >= burn){
      sk_points.col(nEvent-burn) = x%gamma;
      acc_count(nEvent-burn) = acc_prob;
    }
    nEvent++;

    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime";
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        acc_count.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
  }

  List ret ;
  ret["samples"] = sk_points ;
  ret["acc_probs"] = acc_count ;
  return(ret) ;
}

// [[Rcpp::export]]
List hmc_logit_rj_better(double maxTime, const arma::mat& X, const arma::vec& y, arma::vec x0,
                  double epsilon = 0.1, double stoch_time = 1, int nmax = 10^6,
                  int burn = -1, int thin = 1, double prior_sigma2 = 1, double ppi = 0.5,
                  double rw_sig2 = 0.5, double prob_rj = 0.1){

  arma::wall_clock timer;
  timer.tic();

  int p = x0.size(), nEvent = 1;
  int L = std::ceil(stoch_time/epsilon), g_i = 1, move = 1;
  double eps = 1e-10, acc_prob, Ux, Ux_prop, Uz, Uz_prop, prior_diff_gamma;
  double s_gamma, ratio_val, p_gamma_change, q_norm_rj, adj_alpha;

  arma::mat sk_points(p,nmax);
  arma::vec acc_count(nmax), x = x0, x_prop(p), z_vals(p), z_dash(p), gamma(p), gamma_prop(p);
  gamma.zeros();

  if( burn < 0){
    burn = 0;   sk_points.col(0) = x0;
  }
  arma::uvec inds_off_hp = arma::find(arma::abs(x) > eps);
  gamma(inds_off_hp).ones();
  Ux = ll_logit(X, y, x, gamma, prior_sigma2);
  arma::vec x_grad(p);

  arma::uvec inds_on_hp = arma::find(arma::abs(gamma) < eps);

  while( nEvent < nmax + burn){

    for( int it = 0; it < thin; it++){

      // Determine move: 1 = within, 2 = RJ
      move = 1;
      if(R::runif(0,1) < prob_rj){
        move = 2;
      }

      // Within Model
      if(move == 1){
        // Rcout << "move 1";
        if(sum(gamma) > eps){
          z_vals = rnorm(sum(gamma));
          x_prop = x;
          Uz = arma::dot(z_vals,z_vals)/2.0;
          z_dash = z_vals;
          x_grad.zeros();
          for( int l_val = 0; l_val < L; l_val++){
            x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
            z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
            x_prop(inds_off_hp) += epsilon * z_dash;

            x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
            z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
          }
          z_dash = -z_dash;
          Uz_prop = arma::dot(z_dash,z_dash)/2.0;
          Ux_prop = ll_logit(X, y, x_prop, gamma, prior_sigma2);

          acc_prob = std::min(1.0, exp(Ux - Ux_prop + Uz - Uz_prop));
          if( R::runif(0,1) < acc_prob ){
            x = x_prop;
            Ux = Ux_prop;
          }
        }
      } else {
        x_prop = x;
        gamma_prop = gamma;
        g_i = floor(p*R::runif(0,1));
        // Rcout << " move 2: " << g_i << " ";

        // Birth
        if(gamma[g_i] < 0.5){

          gamma_prop[g_i] = 1.0;
          x_prop[g_i] = R::rnorm(0, sqrt(rw_sig2));
          prior_diff_gamma = log(1-ppi) - log(ppi);

          q_norm_rj = -R::dnorm(x_prop[g_i], 0.0, sqrt(rw_sig2), true);

        } else {
          // Death
          gamma_prop[g_i] = 0.0;
          prior_diff_gamma = log(ppi) - log(1-ppi);
          q_norm_rj = R::dnorm(x_prop[g_i], 0.0, sqrt(rw_sig2), true);
        }

        // Rcout << "gam: " << sum(gamma_prop) << " ";

        //Maybe missing a sclaing factor ???
        // ratio_val = -ll_logit(X, y, x_prop, gamma_prop, prior_sigma2)+ll_logit(X, y, x, gamma, prior_sigma2);
        Ux_prop = ll_logit(X, y, x_prop, gamma_prop, prior_sigma2);
        ratio_val = Ux - Ux_prop;
        p_gamma_change = exp(ratio_val-prior_diff_gamma + q_norm_rj);

        if(R::runif(0,1) < p_gamma_change){
          x = x_prop;
          gamma = gamma_prop;
          Ux = Ux_prop;

          inds_off_hp = arma::find(arma::abs(gamma) >= eps);
          inds_on_hp = arma::find(arma::abs(gamma) < eps);
        }
      }
    }

    if( nEvent >= burn){
      sk_points.col(nEvent-burn) = x%gamma;
      acc_count(nEvent-burn) = acc_prob;
    }
    nEvent++;

    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime";
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        acc_count.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
  }

  List ret ;
  ret["samples"] = sk_points ;
  ret["acc_probs"] = acc_count ;
  return(ret) ;
}

// [[Rcpp::export]]
List hmc_logit_rj2(double maxTime, const arma::mat& X, const arma::vec& y, arma::vec x0,
                  double epsilon = 0.1, double stoch_time = 1, int nmax = 10^6,
                  int burn = -1, int thin = 1, double prior_sigma2 = 1,
                  double ppi = 0.5, double rw_sig2 = 0.5, double prob_rj = 0.5){

  arma::wall_clock timer;
  timer.tic();

  int p = x0.size(), nEvent = 1, rj_birth = 1/3, rj_death = 1/3;
  int L = std::ceil(stoch_time/epsilon), g_i = 1, move = 1;
  double eps = 1e-10, acc_prob, Ux, Ux_prop, Uz, Uz_prop, prior_diff_gamma;
  double s_gamma, ratio_val, p_gamma_change, q_norm_rj, adj_alpha;

  arma::mat sk_points(p,nmax);
  arma::vec acc_count(nmax), x = x0, x_prop(p), z_vals(p), z_dash(p), gamma(p), gamma_prop(p);
  gamma.zeros();

  if( burn < 0){
    burn = 0;   sk_points.col(0) = x0;
  }
  arma::uvec inds_off_hp = arma::find(arma::abs(x) > eps);
  gamma(inds_off_hp).ones();
  Ux = ll_logit(X, y, x, gamma, prior_sigma2);
  arma::vec x_grad(p);

  arma::uvec inds_on_hp = arma::find(arma::abs(gamma) < eps);

  while( nEvent < nmax + burn){

    for( int it = 0; it < thin; it++){
      // Determine move: 1 = within, 2 = RJ
      move = 1;
      if(R::runif(0,1) < 0.5){
        move = 2;
      }

      // Within Model
      if(move == 1){
        if(sum(gamma) > eps){
          z_vals = rnorm(sum(gamma));
          x_prop = x;
          Uz = arma::dot(z_vals,z_vals)/2.0; // Momentum Grad (q)
          z_dash = z_vals;
          x_grad.zeros();
          x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;

          // half step momentum first
          z_dash -= epsilon/2.0 * x_grad(inds_off_hp);

          for( int l_val = 0; l_val < L; l_val++){
            // full position move
            x_prop(inds_off_hp) += epsilon * z_dash;
            x_grad(inds_off_hp) = get_grad_logit(X, X.cols(inds_off_hp)*x_prop(inds_off_hp), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
            //full momentum
            if( l_val < L-1){
              z_dash -= epsilon * x_grad(inds_off_hp);
            }
          }
          //Final half step
          z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
          z_dash = -z_dash;
          Uz_prop = arma::dot(z_dash,z_dash)/2.0;
          Ux_prop = ll_logit(X, y, x_prop, gamma, prior_sigma2);

          acc_prob = std::min(1.0, exp(Ux - Ux_prop + Uz - Uz_prop));
          // Rcout << x_prop(0) << " " << x_prop(1)<< "\n";
          if( R::runif(0,1) < acc_prob ){
            x = x_prop;
            Ux = Ux_prop;
          }
          // z_vals = rnorm(sum(gamma));
          // x_prop = x;
          // Uz = arma::dot(z_vals,z_vals)/2.0; // Momentum Grad (q)
          // z_dash = z_vals;
          // x_grad.zeros();
          //
          // for( int l_val = 0; l_val < L; l_val++){
          //
          //   x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
          //   z_dash -= epsilon/2.0 * x_grad(inds_off_hp); // half step off for momentum
          //   x_prop(inds_off_hp) += epsilon * z_dash;     // full position
          //
          //   x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
          //   z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
          // }
          // z_dash = -z_dash;
          // Uz_prop = arma::dot(z_dash,z_dash)/2.0;
          // Ux_prop = ll_logit(X, y, x_prop, gamma, prior_sigma2);
          //
          // acc_prob = std::min(1.0, exp(Ux - Ux_prop + Uz - Uz_prop));
          // Rcout << x_prop(0) << " " << x_prop(1)<< "\n";
          // if( R::runif(0,1) < acc_prob ){
          //   x = x_prop;
          //   Ux = Ux_prop;
          // }
        }
      } else {

        // RJ prop
        // Sample gamma
        for(int j = 0; j < p; j++){
          x_prop = x;
          gamma_prop = gamma;

          if(gamma[j] > 0.0){
            gamma_prop[j] = 0.0;
            prior_diff_gamma = log(ppi) - log(1-ppi);
            q_norm_rj = R::dnorm(x_prop[j], 0.0, sqrt(rw_sig2), true);
          } else {
            gamma_prop[j] = 1.0;
            x_prop[j] = R::rnorm(0, sqrt(rw_sig2));
            prior_diff_gamma = log(1-ppi) - log(ppi);
            q_norm_rj = -R::dnorm(x_prop[inds_on_hp[g_i]], 0.0, sqrt(rw_sig2), true);
          }

          ratio_val = -ll_logit(X, y, x_prop, gamma_prop, prior_sigma2)+ll_logit(X, y, x, gamma, prior_sigma2);
          p_gamma_change = exp(ratio_val-prior_diff_gamma + q_norm_rj);

          if(R::runif(0,1) < p_gamma_change){
            x = x_prop;
            gamma = gamma_prop;
          }
        }

        Ux = ll_logit(X, y, x, gamma, prior_sigma2);
        inds_off_hp = arma::find(arma::abs(gamma) >= eps);
        inds_on_hp = arma::find(arma::abs(gamma) < eps);
      }
    }

    if( nEvent >= burn){
      sk_points.col(nEvent-burn) = x%gamma;
      acc_count(nEvent-burn) = acc_prob;
    }
    nEvent++;

    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime";
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        acc_count.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
  }

  List ret ;
  ret["samples"] = sk_points ;
  ret["acc_probs"] = acc_count ;
  return(ret) ;
}


// [[Rcpp::export]]
List hmc_logit(double maxTime, const arma::mat& X, const arma::vec& y, arma::vec x0,
         double epsilon = 0.1, double stoch_time = 1, int nmax = 10^6,
         int burn = -1, int thin = 1, double prior_sigma2 = 1, double ppi = 0.5){

  arma::wall_clock timer;
  timer.tic();

  int p = x0.size(), nEvent = 1;
  int L = std::ceil(stoch_time/epsilon);
  double eps = 1e-10, acc_prob, Ux, Ux_prop, Uz, Uz_prop, prior_diff_gamma;
  double s_gamma, ratio_val, p_gamma_change;

  arma::mat sk_points(p,nmax);
  arma::vec acc_count(nmax), x = x0, x_prop(p), z_vals(p), z_dash(p), gamma(p), gamma_prop(p);
  gamma.zeros();

  if( burn < 0){
    burn = 0;   sk_points.col(0) = x0;
  }
  arma::uvec inds_off_hp = arma::find(arma::abs(x) > eps);
  gamma(inds_off_hp).ones();
  Ux = ll_logit(X, y, x, gamma, prior_sigma2);
  arma::vec x_grad(p);

  while( nEvent < nmax + burn){

    // Sample gamma
    for(int j = 0; j < p; j++){
      gamma_prop = gamma;
      if(gamma[j] > 0.0){
        gamma_prop[j] = 0.0;
        prior_diff_gamma = -log(ppi) + log(1-ppi);
      } else {
        gamma_prop[j] = 1.0;
        prior_diff_gamma = log(ppi) - log(1-ppi);
      }
      ratio_val = ll_logit(X, y, x, gamma_prop, prior_sigma2) - Ux;
      p_gamma_change = 1.0/(1.0+exp(ratio_val-prior_diff_gamma));

      if(R::runif(0,1) < p_gamma_change){
        gamma = gamma_prop;
        Ux = ratio_val + Ux;
      }
    }

    // Within model moves...
    s_gamma = sum(gamma);
    if( s_gamma > 0 ){
      inds_off_hp = arma::find(arma::abs(gamma) > eps);
      z_vals = rnorm(s_gamma);
      x_prop = x;

      Uz = arma::dot(z_vals,z_vals)/2.0;
      z_dash = z_vals;
      x_grad.zeros();
      for( int l_val = 0; l_val < L; l_val++){
        x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
        z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
        x_prop(inds_off_hp) += epsilon * z_dash;

        x_grad(inds_off_hp) = get_grad_logit(X, X*(x_prop%gamma), y, inds_off_hp) + x(inds_off_hp)/prior_sigma2;
        z_dash -= epsilon/2.0 * x_grad(inds_off_hp);
      }
      z_dash = -z_dash;
      Uz_prop = arma::dot(z_dash,z_dash)/2.0;
      Ux_prop = ll_logit(X, y, x_prop, gamma, prior_sigma2);

      acc_prob = std::min(1.0, exp(Ux - Ux_prop + Uz - Uz_prop));
      if( R::runif(0,1) < acc_prob ){
        x = x_prop;
        Ux = Ux_prop;
      }
    }

    if( nEvent >= burn){
      sk_points.col(nEvent-burn) = x%gamma;
      acc_count(nEvent-burn) = acc_prob;
    }
    nEvent++;

    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime";
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        acc_count.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }
  }

  List ret ;
  ret["samples"] = sk_points ;
  ret["acc_probs"] = acc_count ;
  return(ret) ;
}

#endif
