// --------------------------------------------------------
// gamma2: mean, sd
// --------------------------------------------------------

/**
  * Log gamma density given mean and sd on unit scale
  *
  * @param mean vector of means
  *
  * @param sd standard deviation
  *
  * @return The log of the gamma density of y
  */
real gamma2_lpdf(vector y, vector mean, real sd) {
  int n = num_elements(y);
  vector[n] alpha = (mean / sd)^2;
  vector[n] beta = mean / (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

/**
  * Log gamma density given mean and sd on unit scale
  *
  * @param mean vector of means
  *
  * @param sd vector of standard deviations
  *
  * @return The log of the gamma density of y
  */
real gamma2_lpdf(vector y, vector mean, vector sd) {
  int n = num_elements(y);
  vector[n] alpha = (mean ./ sd)^2;
  vector[n] beta = mean ./ (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

/**
  * Log density of the sum of N i.i.d gamma RVs
  * given mean and sd on unit scale
  *
  * @param mean mean
  *
  * @param sd standard deviation
  *
  * @return The log of the gamma sum density of y
  */
real gamma2_sum_lpdf(vector y, real mean, real sd, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N * ((mean / sd)^2);
  real beta = mean / (sd^2);
  return gamma_lpdf(y | alpha, beta);
}

// --------------------------------------------------------
// gamma3: mean, cv
// --------------------------------------------------------

/**
  * Log gamma density given mean and cv on unit scale
  *
  * @param mean vector of means
  *
  * @param cv coefficient of variation
  *
  * @return The log of the gamma density of y
  */
real gamma3_lpdf(vector y, vector mean, real cv) {
  int n = num_elements(y);
  real alpha = 1 / (cv^2);
  vector[n] beta = 1 / (mean * (cv^2));
  return gamma_lpdf(y | alpha, beta);
}

/**
* Log gamma density given mean and cv on unit scale
*
* @param mean vector of means
*
* @param cv vector of coefficients of variation
*
* @return The log of the gamma density of y
*/
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
