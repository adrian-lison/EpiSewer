// --------------------------------------------------------
// normal: mean, sd
// --------------------------------------------------------

/**
  * Generate truncated normal variate given mean and sd
  * and lower truncation
  *
  * @param mean mean
  *
  * @param sd standard deviation
  *
  * @param lb lower bound
  *
  * @return Truncated normal variate
  */
real normal_lb_rng(real mean, real sd, real lb) {
  real p = normal_cdf(lb | mean, sd);
  real u = uniform_rng(p, 1);
  return mean + sd * inv_Phi(u);
}

/**
  * Generate truncated normal variate given mean and sd
  * and lower truncation
  *
  * @param mean vector of means
  *
  * @param sd standard deviation
  *
  * @param lb lower bound
  *
  * @return Vector of truncated normal variates
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
  * Generate truncated normal variate given mean and sd
  * and lower truncation
  *
  * @param mean vector of means
  *
  * @param sd vector of standard deviations
  *
  * @param lb lower bound
  *
  * @return Vector of truncated normal variates
  */
vector normal_lb_rng(vector mean, vector sd, real lb) {
  int n = num_elements(mean);
  vector[n] p;
  for (i in 1:n) {
    p[i] = normal_cdf(lb | mean[i], sd[i]);
  }
  vector[n] u = to_vector(uniform_rng(p, 1));
  return mean + sd .* inv_Phi(u);
}

/**
  * Generate n standard normal variates
  *
  * @param n Number of i.i.d. samples
  *
  * @return Vector of standard normal variates
  */
vector std_normal_n_rng(int n) {
  return(to_vector(normal_rng(rep_vector(0, n), 1)));
}

// --------------------------------------------------------
// normal2: mean, cv
// --------------------------------------------------------

/**
  * The log of a normal density given mean and coefficient of variation (cv)
  *
  * @param y vector with observed data
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation

  * @return The log of the normal density of y
  */
real normal2_lpdf(vector y, vector mean, vector cv) {
  int n = num_elements(y);
  vector[n] sigma = mean .* cv;
  real tar = normal_lpdf(y | mean, sigma);
  return tar;
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
  real tar = normal_lpdf(y | mean, sigma) - normal_lccdf(lb | mean, sigma);
  for (i in 1:n) {
    if (y[i] < lb) {
      tar += negative_infinity();
      }
  }
  return tar;
}

/**
  * Generate normal variate given mean and coefficient of variation
  *
  * @param mean vector of means
  *
  * @param cv vector of coefficients of variation
  *
  * @return Vector of normal variates with mean exp(mean_log)
  * and standard deviation exp(sd_log)
  */
vector normal2_rng(vector mean, vector cv) {
  int n = num_elements(mean);
  vector[n] sigma = mean .* cv;
  return to_vector(normal_rng(mean, sigma));
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
vector normal2_rng(vector mean, vector cv, real lb) {
  int n = num_elements(mean);
  vector[n] sigma = mean .* cv;
  vector[n] p;
  for (i in 1:n) {
    p[i] = normal_cdf(lb | mean[i], sigma[i]);
  }
  vector[n] u = to_vector(uniform_rng(p, 1));
  return mean + sigma .* inv_Phi(u);
}

// --------------------------------------------------------
// normal prior: mean, sd
// --------------------------------------------------------

/**
  * Normal prior on a parameter, with the option to fix the parameter
  * (i.e. no sampling) by providing a prior with zero variance.
  *
  * @param y Array with the parameter. If the prior has zero variance, the no
  * parameter will be sampled, hence the array has length 0. Otherwise, the
  * array has length 1.
  *
  * @param mean Mean of the prior.
  *
  * @param sd Standard deviation of the prior. If this is zero, then no
  * parameter will be sampled and the mean of the prior will be used instead.
  *
  * @return The log of the prior probability of y
  */
real normal_prior_lpdf(array[] real y, real mean, real sd) {
  if (sd == 0) {
    return(0); // parameter fixed, not sampled
  } else {
    return (normal_lpdf(y | mean, sd));
  }
}

/**
  * Lower truncated normal prior on a parameter, with the option to fix the
  * parameter (i.e. no sampling) by providing a prior with zero variance.
  *
  * @param y Array with the parameter. If the prior has zero variance, the no
  * parameter will be sampled, hence the array has length 0. Otherwise, the
  * array has length 1.
  *
  * @param mean Mean of the prior.
  *
  * @param sd Standard deviation of the prior. If this is zero, then no
  * parameter will be sampled and the mean of the prior will be used instead.
  *
  * @param lb lower bound
  *
  * @return The log of the prior probability of y
  */
real normal_prior_lb_lpdf(array[] real y, real mean, real sd, real lb) {
  if (sd == 0) {
    return(0); // parameter fixed, not sampled
  } else {
    int n = num_elements(y);
    return (normal_lpdf(y | mean, sd) - n * normal_lccdf(lb | mean, sd));
  }
}

/**
  * Upper truncated normal prior on a parameter, with the option to fix the
  * parameter (i.e. no sampling) by providing a prior with zero variance.
  *
  * @param y Array with the parameter. If the prior has zero variance, the no
  * parameter will be sampled, hence the array has length 0. Otherwise, the
  * array has length 1.
  *
  * @param mean Mean of the prior.
  *
  * @param sd Standard deviation of the prior. If this is zero, then no
  * parameter will be sampled and the mean of the prior will be used instead.
  *
  * @param ub upper bound
  *
  * @return The log of the prior probability of y
  */
real normal_prior_ub_lpdf(array[] real y, real mean, real sd, real ub) {
  if (sd == 0) {
    return(0); // parameter fixed, not sampled
  } else {
    int n = num_elements(y);
    return (normal_lpdf(y | mean, sd) - n * normal_lcdf(ub | mean, sd));
  }
}

/**
  * Normal prior on a parameter, with the option to fix the parameter
  * (i.e. no sampling) by providing a prior with zero variance.
  *
  * @param y Array with the parameter. If the prior has zero variance, the no
  * parameter will be sampled, hence the array has length 0. Otherwise, the
  * array has length 1.
  *
  * @param prior The prior for the parameter. This assumes that the prior is
  * stored in an array of length 2, where the first element contains the mean
  * and the second element the standard deviation of the prior.
  *
  * @return The log of the prior probability of y
  */
real normal_prior_lpdf(array[] real y, array[] real prior) {
  return (normal_prior_lpdf(y | prior[1], prior[2]));
}

/**
  * Lower truncated normal prior on a parameter, with the option to fix the
  * parameter (i.e. no sampling) by providing a prior with zero variance.
  *
  * @param y Array with the parameter. If the prior has zero variance, the no
  * parameter will be sampled, hence the array has length 0. Otherwise, the
  * array has length 1.
  *
  * @param prior The prior for the parameter. This assumes that the prior is
  * stored in an array of length 2, where the first element contains the mean
  * and the second element the standard deviation of the prior.
  *
  * @param lb lower bound
  *
  * @return The log of the prior probability of y
  */
real normal_prior_lb_lpdf(array[] real y, array[] real prior, real lb) {
  return (normal_prior_lb_lpdf(y | prior[1], prior[2], lb));
}

/**
  * Upper truncated normal prior on a parameter, with the option to fix the
  * parameter (i.e. no sampling) by providing a prior with zero variance.
  *
  * @param y Array with the parameter. If the prior has zero variance, the no
  * parameter will be sampled, hence the array has length 0. Otherwise, the
  * array has length 1.
  *
  * @param prior The prior for the parameter. This assumes that the prior is
  * stored in an array of length 2, where the first element contains the mean
  * and the second element the standard deviation of the prior.
  *
  * @param ub upper bound
  *
  * @return The log of the prior probability of y
  */
real normal_prior_ub_lpdf(array[] real y, array[] real prior, real ub) {
  return (normal_prior_ub_lpdf(y | prior[1], prior[2], ub));
}

/**
  * Calculate the mean of a truncated normal distribution (truncated below zero)
  * with mean mu and standard deviation sigma.
  *
  * @param mu Mean of the untruncated normal distribution
  *
  * @param sigma Standard deviation of the untruncated normal distribution
  *
  * @return The mean of the truncated normal distribution
  */
real trunc_normal_mean(real mu, real sigma) {
  if (sigma == 0) {
    return mu;
  } else {
    real alpha = -mu / sigma;
    real phi_alpha = exp(std_normal_lpdf(alpha));
    real Phi_alpha = Phi_approx(alpha);
    return mu + sigma * phi_alpha / (1 - Phi_alpha);
  }
}
