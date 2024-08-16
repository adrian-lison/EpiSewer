/**
* Density of a log-gamma random variable (L. J. Halliwell 2018)
*/
real loggamma_lpdf(vector y, real alpha, real beta){
    return(sum(-lgamma(alpha) - (exp(y) * beta) + (y + log(beta)) * alpha));
}

real loggamma_lpdf(vector y, real alpha, vector beta){
    return(sum(-lgamma(alpha) - (exp(y) .* beta) + (y + log(beta)) * alpha));
}

real loggamma_lpdf(vector y, vector alpha, real beta){
    return(sum(-lgamma(alpha) - (exp(y) * beta) + (y + log(beta)) .* alpha));
}

real loggamma_lpdf(vector y, vector alpha, vector beta){
    return(sum(-lgamma(alpha) - (exp(y) .* beta) + (y + log(beta)) .* alpha));
}

// Mean and CV parameterization
real loggamma3_lpdf(vector y, real mean, real cv){
  real alpha = 1 / (cv^2);
  real beta = 1 / (mean * (cv^2));
  return loggamma_lpdf(y | alpha, beta);
}

real loggamma3_lpdf(vector y, vector mean, real cv){
  int n = num_elements(y);
  real alpha = 1 / (cv^2);
  vector[n] beta = 1 / (mean * (cv^2));
  return loggamma_lpdf(y | alpha, beta);
}
