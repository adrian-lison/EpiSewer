/**
  * Discretize a continuous probability distribution on the log scale
  *
  * Discretizes a continous probability distribution to integer buckets on the
  * log scale under the assumption of single censoring.
  *
  * This code was adapted from
  * https://github.com/epiforecasts/EpiNow2
  * and is itself based on code from
  * https://github.com/epiforecasts/epinowcast
  * Copyright: EpiNow2 and epinowcast authors
  *
  * @param params The parameters of the distribution
  *
  * @param max_x The maximum value to consider for discretization
  *
  * @param dist_type The type of distribution to use:
  *   1 for exponential, 2 for gamma, 3 for lognormal
  *
  * @return A vector of log probabilities for each value from 0 to max_x
  */
vector discretise_dist_log(int dist_type, real mean, real cv, int max_x) {
  int n = max_x + 1;
  vector[n] lpmf;
  vector[n] upper_lcdf;
  if (dist_type == 0) {
    reject(
      "Dist_type = 0 corresponds to an already discretized ",
      "non-parametric function and is not supported."
      );
  } else if (dist_type == 1) {
    real lambda = 1 / mean;
    for (i in 1:n) {
      upper_lcdf[i] = exponential_lcdf(i | lambda);
    }
  } else if (dist_type == 2) {
    real alpha = (1/cv)^2;
    real beta = alpha / mean;
    for (i in 1:n) {
      upper_lcdf[i] = gamma_lcdf(i | alpha, beta);
    }
  } else if (dist_type == 3) {
    real sigma2 = log(cv^2 + 1);
    real mu = log(mean) - sigma2 / 2;
    for (i in 1:n) {
      upper_lcdf[i] = lognormal_lcdf(i | mu, sqrt(sigma2));
    }
  } else {
    reject("Unknown distribution function provided.");
  }
  // discretise
  if (n > 1) {
    lpmf[1] = upper_lcdf[1];
    lpmf[2:n] = log_diff_exp(upper_lcdf[2:n], upper_lcdf[1:(n-1)]);
    // normalize
    lpmf = lpmf - upper_lcdf[n];
  } else {
    lpmf[1] = 0;
  }
  return(lpmf);
}

/**
  * Discretize a continuous probability distribution
  *
  * Discretizes a continous probability distribution to integer buckets under
  * the assumption of single censoring.
  *
  * This code was adapted from
  * https://github.com/epiforecasts/EpiNow2
  * and is itself based on code from
  * https://github.com/epiforecasts/epinowcast
  * Copyright: EpiNow2 and epinowcast authors
  *
  * @param params The parameters of the distribution
  *
  * @param max_x The maximum value to consider for discretization
  *
  * @param dist_type The type of distribution to use:
  *   1 for exponential, 2 for gamma, 3 for lognormal
  *
  * @return A vector of probabilities for each value from 0 to max_x
  */
vector discretise_dist(int dist_type, real mean, real cv, int max_x) {
  return(exp(discretise_dist_log(dist_type, mean, cv, max_x)));
}
