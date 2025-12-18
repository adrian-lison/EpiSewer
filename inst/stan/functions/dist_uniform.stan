// --------------------------------------------------------
// uniform: a, b
// --------------------------------------------------------

real uniform_prior_lpdf(array[] real y, real a, real b) {
  if (a == b) {
    return(0); // parameter fixed, not sampled
  } else {
    return uniform_lpdf(y | a, b);
  }
}

real uniform_prior_lpdf(array[] real y, array[] real prior) {
  return uniform_prior_lpdf(y | prior[1], prior[2]);
}
