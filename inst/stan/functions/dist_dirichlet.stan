/**
* Transform unconstrained vector to simplex
*/
vector transform_to_simplex(vector z) {
  vector[num_elements(z)] g = exp(z);
  vector[num_elements(z)] theta = g / sum(g);
  return theta;
}

real dirichlet_raw_lpdf(vector z, vector alpha) {
  real lpdf = 0;
  lpdf += gamma_lpdf(exp(z) | alpha, 1.0);  // Gamma(alpha, 1)
  lpdf += sum(z);  // Jacobian adjustment for exp()
  return(lpdf);
}
