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

real gamma3_sum_lpdf(vector y, real mean, real cv, vector N) {
  int n = num_elements(y);
  vector[n] alpha = N / (cv^2); // sum of gammas with same shape and scale
  real beta = 1 / (mean * (cv^2));
  return gamma_lpdf(y | alpha, beta);
}
