// Sum of gamma distributed RVs
real gamma3_sum_lpdf(vector y, real mean, real cv, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N / (cv^2); // sum of gammas with same shape and scale
  real beta = 1 / (mean * (cv^2));
  return gamma_lpdf(y | alpha, beta);
}

// Log of sum of gamma distributed RVs
real gamma3_sum_log_lpdf(vector y, real mean, real cv, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N / (cv^2); // sum of gammas with same shape and scale
  real beta = 1 / (mean * (cv^2));
  return loggamma_lpdf(y | alpha, beta);
}

// Non-centered paramaterization of a normal approximation for the
// sum of N i.i.d. Gamma distributed RVs with mean 1 and a specified cv
vector gamma_sum_approx(real cv, vector N, vector noise_noncentered) {
  // sqrt(N) * cv is the standard deviation of the sum of Gamma distributions
  return N + noise_noncentered .* sqrt(N) * cv;
}

// Non-centered paramaterization of a normal approximation for the
// log of the sum of N i.i.d. Gamma distributed RVs with mean 1 and a specified cv
vector gamma_sum_log_approx(real cv, vector N, vector noise_noncentered) {
  // sqrt(N) * cv is the standard deviation of the sum of Gamma distributions
  return log(softplus(N + noise_noncentered .* sqrt(N) * cv, 10));
}
