
/**
  * Beta prior on a parameter, parameterized via mean and sd, with the option
  * to fix the parameter (i.e. no sampling) by providing a prior with zero
  * variance.
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
real beta2_prior_lpdf(array[] real y, real mean, real sd) {
  if (sd == 0) {
    return(0); // parameter fixed, not sampled
  } else {
    real alpha = (mean*(1-mean)/(sd^2)-1)*mean;
    real beta = alpha * (1-mean)/mean;
    return (beta_lpdf(y | alpha, beta));
  }
}

/**
  * Beta prior on a parameter, parameterized via mean and sd, with the option
  * to fix the parameter (i.e. no sampling) by providing a prior with zero
  * variance.
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
real beta2_prior_lpdf(array[] real y, array[] real prior) {
  return (beta2_prior_lpdf(y | prior[1], prior[2]));
}
