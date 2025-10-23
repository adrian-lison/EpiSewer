/**
 * Utility functions for Hilbert space approximate Bayesian Gaussian processes
 * as described in Riutort-Mayol et al. 2023.
 * This implementation is strongly inspired by the implementation in the
 * EpiNow2 package, credits to the EpiNow2 authors.
 *
 * - Riutort-Mayol et al. 2023 paper: https://doi.org/10.1007/s11222-022-10167-2
 * - EpiNow2 package: https://github.com/epiforecasts/EpiNow2/
 */

/**
  * Compute spectral densities for Matern kernel
  *
  * @param nu Smoothness
  * @param sigma Magnitude
  * @param rho Length scale
  * @param L Boundary condition
  * @param m Number of basis functions
  * @return A vector of spectral densities
  */
vector diagSPD_Matern(real nu, real sigma, real rho, real L, int m) {
  real numerator;
  vector[m] denominator;
  // vector of frequency locations
  vector[m] omega = pi() / (2 * L) * linspaced_vector(m, 1, m);
  if (nu == 0.5) {
    // Matern 1/2 kernel
    numerator = 2;
    denominator = rho * ((1 / rho)^2 + omega^2);
  } else if (nu == 1.5) {
    // Matern 3/2 kernel
    numerator = 4 * (sqrt(3) / rho)^3;
    denominator = pow((sqrt(3) / rho)^2 + omega^2, 2);
  } else if (nu == 2.5) {
    // Matern 5/2 kernel
    numerator = 3 * (sqrt(5) / rho)^5;
    denominator = 2 * pow((sqrt(5) / rho)^2 + omega^2, 3);
  } else {
    reject("nu must be one of 1/2, 3/2 or 5/2; found nu=", nu);
  }
  return sigma * sqrt(numerator ./ denominator);
}

/**
  * Compute eigenvalue matrix of Laplacian operator for Hilbert space approx
  *
  * @param n Number of points to be evaluated
  * @param m Number of basis functions
  * @param L Boundary condition
  * @param x Points in input space
  *
  * @return A matrix of eigenvalues
  */
matrix gp_eigenvalue_matrix(int n, int m, real L, vector x) {
  matrix [n, m] PHI = inv_sqrt(L) * sin(
    diag_post_multiply(
      rep_matrix(pi() / (2 * L) * (x + L), m),
      linspaced_vector(m, 1, m)
    ));
  return PHI;
}
