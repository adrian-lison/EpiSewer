/**
  * The exponentially modified gaussian density given mean, sd and shape parameters
  *
  * @param y vector with observed data
  * 
  * @param mean vector of means
  *
  * @param sd standard deviation (same for all entries)
  *
  * @param shape shape (same for all entries)
  * 
  * @return The exponentially modified gaussian density of y
  */
real exp_mod_normal2_lpdf(vector y, vector mean, real sd, real shape) {
  int n = num_elements(y);
  real sigma = sd/sqrt(1+shape^2);
  real lambda = 1/(shape*sigma);
  vector[n] mu = mean - shape * sigma;
  return exp_mod_normal_lpdf(y | mu, sigma, lambda);
}

/**
  * Generate exponentially modified gaussian variate given mean, sd and shape parameters
  * 
  * @param y vector with observed data
  * 
  * @param mean vector of means
  *
  * @param sd standard deviation (same for all entries)
  *
  * @param shape shape (same for all entries)
  * 
  * @return Vector of exponentially modified gaussian variates
  */
array[] real exp_mod_normal2_rng(vector mean, real sd, real shape) {
  int n = num_elements(y);
  real sigma = sd/sqrt(1+shape^2);
  real lambda = 1/(shape*sigma);
  vector[n] mu = mean - shape * sigma;
  return exp_mod_normal_rng(mu, sigma, lambda);
}