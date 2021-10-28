/* 
   Bayesian version of the TVGEV function

   For now, missing values are not allowed.

*/

functions{

  real gev_lpdf(vector y, vector mu, vector sigma, vector xi) {
    vector[rows(y)] t;
    vector[rows(y)] lp;
    int n;
    n = rows(y);
    for(i in 1:n){
      t[i] = xi[i]==0 ? exp((mu[i] - y[i]) / sigma[i]) :
	pow(1.0 + xi[i] * ((y[i] - mu[i]) / sigma[i]), -1.0 / xi[i]);
      lp[i] = -log(sigma[i]) + (xi[i] + 1.0) * log(t[i]) - t[i];
    }
    return sum(lp);
  }
  
  real gev_lcdf(vector y, vector mu, vector sigma, vector xi) {
    vector[rows(y)] t;
    int n;
    n = rows(y);
    for(i in 1:n){
      t[i] = xi[i]==0 ? exp((mu[i] - y[i]) / sigma[i]) :
	pow(1.0 + xi[i] * ((y[i] - mu[i]) / sigma[i]), -1.0 / xi[i]);
    }
    return -sum(t);
  }
  
  real gev_cdf(real y, real mu, real sigma, real xi) {
    real t;
    t = xi == 0 ? exp((mu - y) / sigma) :
      pow(1.0 + xi * ((y - mu) / sigma), -1.0 / xi);
    return exp(-t);
  }
}
data {
  
  // observations
  int<lower=0> n_obs;
  int<lower=0> n_miss;
  int<lower=0> n_cens;
  
  vector[n_obs] y_obs;
  vector[n_cens] yL_cens;
  vector[n_cens] yU_cens;
  vector[n_cens] code_cens;
  
  // Maybe use an array later?
  int<lower=1> p_mu;
  int<lower=1> p_sigma;
  int<lower=1> p_xi;
  
  // Flags
  int<lower=0> cst_mu;
  int<lower=0> cst_sigma;
  int<lower=0> cst_xi;
  
  // design matrices
  matrix[n_obs, p_mu] X_mu_obs;
  matrix[n_obs, p_sigma] X_sigma_obs;
  matrix[n_obs, p_xi] X_xi_obs;

 
  matrix[n_miss, p_mu] X_mu_miss;
  matrix[n_miss, p_sigma] X_sigma_miss;
  matrix[n_miss, p_xi] X_xi_miss;
  
    
  matrix[n_cens, p_mu] X_mu_cens;
  matrix[n_cens, p_sigma] X_sigma_cens;
  matrix[n_cens, p_xi] X_xi_cens;  
  
  // prior covariances
  matrix[p_mu, p_mu] cov_psi_mu;
  matrix[p_sigma, p_sigma] cov_psi_sigma;
  matrix[p_xi, p_xi] cov_psi_xi;
  
  // prior means
  vector[p_mu] mean_psi_mu;
  vector[p_sigma] mean_psi_sigma;
  vector[p_xi] mean_psi_xi;
  
}

/* Note that we must not sample the censored observations because
   these are coped with by adding their contribution in the
   log-likelihood */

parameters {
  
  vector[p_mu] psi_mu;
  vector[p_sigma] psi_sigma;
  vector[p_xi] psi_xi;
  vector[n_miss] y_miss;
  // vector[n_cens] y_cens;
  
}

model {

  vector[n_obs] mu_obs;
  vector[n_obs] sigma_obs;
  vector[n_obs] xi_obs;

  vector[n_miss] mu_miss;
  vector[n_miss] sigma_miss;
  vector[n_miss] xi_miss;
    
  vector[n_cens] mu_cens;
  vector[n_cens] sigma_cens;
  vector[n_cens] xi_cens;

  // sample from the (Gaussian) priors
  // add a test for a constant parameter 
  psi_mu ~ multi_normal(mean_psi_mu, cov_psi_mu);
  psi_sigma ~ multi_normal(mean_psi_sigma, cov_psi_sigma);
  psi_xi ~ multi_normal(mean_psi_xi, cov_psi_xi);
  
  if (cst_mu) {
    for (i in 1:n_obs) mu_obs[i] = psi_mu[1];
    for (i in 1:n_miss) mu_miss[i] = psi_mu[1];
    for (i in 1:n_cens) mu_cens[i] = psi_mu[1];
  } else {
    mu_obs = X_mu_obs * psi_mu;
    if (n_miss) mu_miss = X_mu_miss * psi_mu;
    if (n_cens) mu_cens = X_mu_cens * psi_mu;
  }
  
  if (cst_sigma) {
    for (i in 1:n_obs) sigma_obs[i] = psi_sigma[1];
    for (i in 1:n_miss) sigma_miss[i] = psi_sigma[1];
    for (i in 1:n_cens) sigma_cens[i] = psi_sigma[1];
  } else {
    sigma_obs = X_sigma_obs * psi_sigma;
    if (n_miss) sigma_miss = X_sigma_miss * psi_sigma;
    if (n_cens) sigma_cens = X_sigma_cens * psi_mu;
  }
  
  if (cst_xi) {
    for (i in 1:n_obs) xi_obs[i] = psi_xi[1];
    for (i in 1:n_miss) xi_miss[i] = psi_xi[1];
    for (i in 1:n_cens) xi_cens[i] = psi_xi[1];
  } else {
    xi_obs = X_xi_obs * psi_xi;
    if (n_miss) xi_miss = X_xi_miss * psi_xi;
    if (n_cens) xi_cens = X_xi_cens * psi_xi;
  }
  
  // likelihood
  y_obs ~ gev(mu_obs, sigma_obs, xi_obs);
  if (n_miss) y_miss ~ gev(mu_miss, sigma_miss, xi_miss);

  // likelihood contributions of the censored observations
  for (i in 1:n_cens) {
    if (code_cens[i] == 1) {
      target += log(1.0 - gev_cdf(yL_cens[i], mu_cens[i], sigma_cens[i], xi_cens[i]));
    } else if (code_cens[i] == 2) {
      target += log(gev_cdf(yU_cens[i], mu_cens[i], sigma_cens[i], xi_cens[i]) -
  		    gev_cdf(yL_cens[i], mu_cens[i], sigma_cens[i], xi_cens[i]));
    } else if (code_cens[i] == 3) {
      target += log(gev_cdf(yU_cens[i], mu_cens[i], sigma_cens[i], xi_cens[i]));
    }
  }
  
}
