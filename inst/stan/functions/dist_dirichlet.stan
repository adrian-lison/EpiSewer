/**
* Parameterization of the Dirichlet distribution via independent Gamma distributions.
*/
real dirichlet_raw_lpdf(vector z, vector alpha) {
  real lpdf = 0;
  lpdf += gamma_lpdf(exp(z) | alpha, 1.0);  // Gamma(alpha, 1)
  lpdf += sum(z);  // Jacobian adjustment for exp()
  return(lpdf);
}

/**
* Parameterization of the Dirichlet distribution via independent Gamma distributions,
* with dimensionality reduction (only n-1 Gammas) to avoid overparameterization
*/
real dirichlet_reduced_lpdf(vector z_raw, vector alpha) {
    int   n = num_elements(alpha);
    vector[n] z = append_row(z_raw, 0);
    vector[n] g = exp(z);
    real  T = sum(g);
    real alpha0 = sum(alpha);
    real log_dir = dot_product(alpha - 1, z) - (alpha0 - n) * log(T); // Dirichlet density
    real log_jac = sum(z_raw) - n * log(T); // Jacobian adjustment for softmax
    return log_dir + log_jac;
  }
