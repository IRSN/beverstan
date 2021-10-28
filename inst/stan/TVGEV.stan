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
  
}
data {
  
  // observations
  int<lower=0> n;
  vector[n] y;
  
  // Maybe use an array later?
  int<lower=1> p_mu;
  int<lower=1> p_sigma;
  int<lower=1> p_xi;
  
  // Flags
  int<lower=0> cst_mu;
  int<lower=0> cst_sigma;
  int<lower=0> cst_xi;
  
  // design matrices
  matrix[n, p_mu] X_mu;
  matrix[n, p_sigma] X_sigma;
  matrix[n, p_xi] X_xi;
  
  // prior covariances
  matrix[p_mu, p_mu] cov_psi_mu0;
  matrix[p_sigma, p_sigma] cov_psi_sigma0;
  matrix[p_xi, p_xi] cov_psi_xi0;
  
  // prior means
  vector[p_mu] mean_psi_mu0;
  vector[p_sigma] mean_psi_sigma0;
  vector[p_xi] mean_psi_xi0;
  
}

parameters {
  
  vector[p_mu] psi_mu;
  vector[p_sigma] psi_sigma;
  vector[p_xi] psi_xi;

}

model {

  vector[n] mu;
  vector[n] sigma;
  vector[n] xi;
  
  // sample from the (Gaussian) priors
  // add a test for a constant parameter 
  psi_mu ~ multi_normal(mean_psi_mu0, cov_psi_mu0);
  if (cst_mu) {
    for (i in 1:n) {
      mu[i] = psi_mu[1];
    }
  } else {
    mu = X_mu * psi_mu;
  }
  
  psi_sigma ~ multi_normal(mean_psi_sigma0, cov_psi_sigma0);
  if (cst_sigma) {
    for (i in 1:n) {
      sigma[i] = psi_sigma[1];
    }
  } else {
    sigma = X_sigma * psi_sigma;
  }
  
  psi_xi ~ multi_normal(mean_psi_xi0, cov_psi_xi0);
  if (cst_xi) {
    for (i in 1:n) {
      xi[i] = psi_xi[1];
    }
  } else {
    xi = X_xi * psi_xi;
  }
  
  // likelihood
  y ~ gev(mu, sigma, xi);
  
}
