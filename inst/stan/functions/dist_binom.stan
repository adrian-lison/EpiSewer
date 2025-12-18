// --------------------------------------------------------
// binom_approx: size n, probability p
// --------------------------------------------------------

/**
  * The log of the binomial density given size and probability, with an
  approximation for approximate
  *
  * @param y observed data
  *
  * @param n size
  *
  * @param p probability
  *
  * @return The log of the continuous binomial density of y
  */
real binom_approx_lpdf(real y, real n, real p) {
  real logC = lgamma(n + 1) - (lgamma(y + 1) + lgamma(n - y + 1));
  real logP = y * log(p) + (n - y) * log(1 - p);
  return(logC + logP);
  }

real binom_approx_lpdf(vector y, vector n, vector p) {
  int N = num_elements(y);
  vector[N] logC = lgamma(n + 1) - (lgamma(y + 1) + lgamma(n - y + 1));
  vector[N] logP = y .* log(p) + (n - y) .* log(1 - p);
  return(sum(logC + logP));
  }

vector binom_approx_lpdfs(vector y, real n, real p) {
  int N = num_elements(y);
  vector[N] logC = lgamma(n + 1) - (lgamma(y + 1) + lgamma(n - y + 1));
  vector[N] logP = y * log(p) + (n - y) * log(1 - p);
  return(logC + logP);
  }

vector binom_approx_lpdfs(vector y, vector n, vector p) {
  int N = num_elements(y);
  vector[N] logC = lgamma(n + 1) - (lgamma(y + 1) + lgamma(n - y + 1));
  vector[N] logP = y .* log(p) + (n - y) .* log(1 - p);
  return(logC + logP);
  }
