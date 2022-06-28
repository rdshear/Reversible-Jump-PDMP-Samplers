// Adapted from https://projecteuclid.org/download/pdfview_1/euclid.ejs/1513306866
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;
  matrix[n,p] x;
  real<lower=1> nu_local ; // degrees of freedom for the half - t priors for lambdas (=1 horseshoe)
  real<lower=0> c ; // slab scale
  // Taus section
  real<lower=0> tau ;
  // real<lower=0> scale_global ; // scale for the half - t prior for tau
  // real<lower=1> nu_global ; // degrees of freedom for the half - t priors for tau
}

parameters {
  // auxiliary variables that define the global and local parameters
  vector[p] aux_z;
  vector<lower=0>[p] lambda ; // local shrinkage
  // real<lower=0> tau ; // global shrinkage
}

transformed parameters {
  vector<lower=0>[p] lambda_tilde ; // truncated local shrinkage
  vector[p] beta ; // regression coefficients
  vector[n] mu; // latent values
  lambda_tilde = sqrt(c^2 * square(lambda) ./(c^2 + tau^2*square(lambda)));
  beta = aux_z .* lambda_tilde * tau ;
  mu = x*beta ;
}

model {
  // half - t priors for lambdas and tau
  aux_z ~ normal(0, 1);
  lambda ~ student_t(nu_local, 0, 1);
  // tau ~ student_t(nu_global, 0, scale_global); // default (0, 1)?

  // robust regression likelihood
  for (i in 1:n){
    target += log_mix(0.5,normal_lpdf(y[i] | mu[i], 1),
    normal_lpdf(y[i] | mu[i], 10));
  }
}
