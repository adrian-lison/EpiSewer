real gamma2_lpdf(vector y, vector mean, real sd) {
  int n = num_elements(y);
  vector[n] alpha = (mean / sd)^2;
  vector[n] beta = mean / (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

real gamma2_sum_lpdf(vector y, real mean, real sd, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N * ((mean / sd)^2);
  real beta = mean / (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

real gamma3_lpdf(vector y, vector mean, real cv) {
  int n = num_elements(y);
  real alpha = 1 / (cv^2);
  vector[n] beta = 1 / (mean * (cv^2));
  return gamma_lpdf(y | alpha, beta);
}

real gamma3_lpdf(vector y, vector mean, vector cv) {
  int n = num_elements(y);
  vector[n] cv_squared = (cv^2);
  vector[n] alpha = 1 / cv_squared;
  vector[n] beta = 1 / (mean .* cv_squared);
  return gamma_lpdf(y | alpha, beta);
}

/**
  * Generate Gamma variate given mean and coefficient of variation
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @return Vector of Gamma variates
  */
vector gamma3_rng(vector mean, vector cv) {
  int n = num_elements(mean);
  vector[n] cv_squared = (cv^2);
  vector[n] alpha = 1 / cv_squared;
  vector[n] beta = 1 / (mean .* cv_squared);
  return to_vector(gamma_rng(alpha, beta));
}

real gamma3_sum_lpdf(vector y, real mean, real cv, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N / (cv^2); // sum of gammas with same shape and scale
  real beta = 1 / (mean * (cv^2));
  return gamma_lpdf(y | alpha, beta);
}

// Non-centered paramaterization of a normal approximation for the
// sum of N i.i.d. Gamma distributed RVs with mean 1 and a specified cv
vector gamma_sum_approx(real cv, vector N, vector noise_noncentered) {
  // sqrt(N) * cv is the standard deviation of the sum of Gamma distributions
  return N + noise_noncentered .* sqrt(N) * cv;
}
