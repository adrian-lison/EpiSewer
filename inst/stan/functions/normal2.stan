/**
  * The log of the normal density given mean and coefficient of variation (cv)
  *
  * @param y vector with observed data
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @return The log of the normal density of y
  */
real normal2_lpdf(vector y, vector mean, vector cv) {
  int n = num_elements(y);
  vector[n] sigma = mean .* cv;
  return normal_lpdf(y | mean, sigma);
}

/**
  * Generate truncated normal variate given mean and coefficient of variation
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @return Vector of normal variates with mean exp(mean_log)
  * and standard deviation exp(sd_log)
  */
vector normal2_lb_rng(vector mean, vector cv, real lb) {
  int n = num_elements(mean);
  vector[n] sigma = mean .* cv;
  vector[n] p;
  for (i in 1:n) {
    p[i] = normal_cdf(lb | mean[i], sigma[i]);
  }
  vector[n] u = to_vector(uniform_rng(p, 1));
  return mean + sigma .* inv_Phi(u);
}

