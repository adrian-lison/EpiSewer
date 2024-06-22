/**
  * Generate truncated normal variate given mean and coefficient of variation
  * and lower truncation
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @param lb lower bound
  *
  * @return Vector of normal variates with mean exp(mean_log)
  * and standard deviation exp(sd_log)
  */
vector normal_lb_rng(vector mean, real sd, real lb) {
  int n = num_elements(mean);
  vector[n] p;
  for (i in 1:n) {
    p[i] = normal_cdf(lb | mean[i], sd);
  }
  vector[n] u = to_vector(uniform_rng(p, 1));
  return mean + sd * inv_Phi(u);
}

/**
  * The log of a normal density given mean and coefficient of variation (cv)
  * with lower truncation
  *
  * @param y vector with observed data
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @param lb lower bound
  *
  * @return The log of the normal density of y
  */
real normal2_lpdf(vector y, vector mean, vector cv, real lb) {
  int n = num_elements(y);
  vector[n] sigma = mean .* cv;
  real tar = normal_lpdf(y | mean, sigma);
  for (i in 1:n) {
    if (y[i] < lb) {
      tar += negative_infinity();
      }
  }
  tar += -normal_lccdf(lb | mean, sigma);
  return tar;
}

/**
  * Generate truncated normal variate given mean and coefficient of variation
  * and lower truncation
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @param lb lower bound
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

vector std_normal_n_rng(int n) {
  return(to_vector(normal_rng(rep_vector(0, n), 1)));
}
