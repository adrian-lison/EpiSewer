/**
  * The log of the lognormal density given mean and sd parameters on log scale
  *
  * @param y vector with observed data
  *
  * @param mean_log vector of log of means
  *
  * @param sd_log log of standard deviation (same for all entries)
  *
  * @return The log of the lognormal density of y with mean exp(mean_log)
  * and standard deviation exp(sd_log)
  */
real lognormal2_lpdf(vector y, vector mean_log, real sd_log) {
  int n = num_elements(y);
  vector[n] sigma2 = log1p_exp(2*(sd_log - mean_log));
  vector[n] mu = mean_log - sigma2/2;
  return lognormal_lpdf(y | mu, sqrt(sigma2));
}

/**
  * Generate lognormal variate given mean and sd parameters on log scale
  *
  * @param mean_log vector of log of means
  *
  * @param sd_log log of standard deviation (same for all entries)
  *
  * @return Vector of lognormal variates with mean exp(mean_log)
  * and standard deviation exp(sd_log)
  */
array[] real lognormal2_rng(vector mean_log, real sd_log) {
  int n = num_elements(mean_log);
  vector[n] sigma2 = log1p_exp(2*(sd_log - mean_log));
  vector[n] mu = mean_log - sigma2/2;
  return lognormal_rng(mu, sqrt(sigma2));
}

/**
  * The log of the lognormal density given mean on log scale and sigma on unit
  * scale
  *
  * @param y vector with observed data
  *
  * @param mean_log vector of log of means
  *
  * @param sigma scale of the lognormal - not the standard deviation!
  *
  * @return The log of the lognormal density of y
  */
real lognormal3_lpdf(vector y, vector mean_log, real sigma) {
  int n = num_elements(y);
  vector[n] mu = mean_log - square(sigma)/2;
  return lognormal_lpdf(y | mu, sigma);
}

/**
  * Generate lognormal variate given mean on log scale and sigma on unit scale
  *
  * @param mean_log vector of log of means
  *
  * @param sigma scale of the lognormal - not the standard deviation!
  *
  * @return Vector of lognormal variates
  */
array[] real lognormal3_rng(vector mean_log, real sigma) {
  int n = num_elements(mean_log);
  vector[n] mu = mean_log - square(sigma)/2;
  return lognormal_rng(mu, sigma);
}

/**
  * The log of the lognormal density given mean on log scale and coefficient of
  * variation on unit scale
  *
  * @param y vector with observed data
  *
  * @param mean_log vector of log of means
  *
  * @param cv coefficient of variation
  *
  * @return The log of the lognormal density of y
  */
real lognormal4_lpdf(vector y, vector mean_log, real cv) {
  int n = num_elements(y);
  real sigma2 = log1p(cv^2);
  vector[n] mu = mean_log - sigma2/2;
  return lognormal_lpdf(y | mu, sqrt(sigma2));
}

/**
  * Generate lognormal variate given mean on log scale and coefficient of
  * variation on unit scale
  *
  * @param mean_log vector of log of means
  *
  * @param cv coefficient of variation
  *
  * @return Vector of lognormal variates
  */
array[] real lognormal4_rng(vector mean_log, real cv) {
  int n = num_elements(mean_log);
  real sigma2 = log1p(cv^2);
  vector[n] mu = mean_log - sigma2/2;
  return lognormal_rng(mu, sqrt(sigma2));
}

/**
  * The log of the density of a lognormal distribution on the log scale (i.e.
  * a normal) with mean on log scale and coefficient of variation on unit scale
  *
  * @param y vector with observed data (on log scale)
  *
  * @param mean_log vector of log of means
  *
  * @param cv coefficient of variation
  *
  * @return The log of the normal density of y
  */
real lognormal_log_lpdf(vector y, vector mean_log, real cv) {
  int n = num_elements(y);
  real sigma2 = log1p(cv^2);
  vector[n] mu = mean_log - sigma2/2;
  return normal_lpdf(y | mu, sqrt(sigma2));
}

/**
  * Generate normal variates corresponding to a lognormal distribution on the
  * log scale, with mean on log scale and coefficient of variation on unit scale
  *
  *
  * @param mean_log vector of log of means
  *
  * @param cv coefficient of variation
  *
  * @return vector with normal variates
  */
array[] real lognormal_log_rng(vector mean_log, real cv) {
  int n = num_elements(mean_log);
  real sigma2 = log1p(cv^2);
  vector[n] mu = mean_log - sigma2/2;
  return normal_rng(mu, sqrt(sigma2));
}
