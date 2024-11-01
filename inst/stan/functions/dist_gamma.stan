// Mean and sd parameterization
real gamma2_lpdf(vector y, vector mean, real sd) {
  int n = num_elements(y);
  vector[n] alpha = (mean / sd)^2;
  vector[n] beta = mean / (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

real gamma2_lpdf(vector y, vector mean, vector sd) {
  int n = num_elements(y);
  vector[n] alpha = (mean ./ sd)^2;
  vector[n] beta = mean ./ (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

real gamma2_sum_lpdf(vector y, real mean, real sd, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N * ((mean / sd)^2);
  real beta = mean / (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

// Mean and CV parameterization
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
  if (min(beta) == 0) {
    reject(
        "Inverse scale parameter of Gamma distribution is zero. | ",
        "Mean: ", mean, " | ",
        "CV: ", cv, " | "
        );
  }
  return to_vector(gamma_rng(alpha, beta));
}
