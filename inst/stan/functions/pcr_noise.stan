/**
  * Compute coefficient of variation (CV) in dPCR quantification as a function
  * of the concentration (no pre-PCR noise)
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV (all other noise factors)
  *
  * @param m total number of partitions (i.e. sum of partitions of all replicates)
  *
  * @param c conversion factor c = sv, i.e. should be partition volume v
  * multiplied by scaling factor s (concentration in assay / concentration in wastewater)
  *
  * @return A vector with the corresponding coefficients of variation
  */
vector cv_dPCR(vector lambda, vector m, real c) {
  int N = num_elements(lambda);
  vector[N] total_var = (exp(lambda * c) - 1) ./ (m * c^2);
  vector[N] cv = sqrt(total_var) ./ lambda;
  return cv;
}


/**
  * Compute coefficient of variation (CV) in dPCR quantification as a function
  * of the concentration, under pre-PCR noise
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV (all other noise factors)
  *
  * @param m total number of partitions (i.e. sum of partitions of all replicates)
  *
  * @param c conversion factor c = sv, i.e. should be partition volume v
  * multiplied by scaling factor s (concentration in assay / concentration in original sample)
  *
  * @param pre_type type of pre-PCR noise (0: gamma, 1: log-normal)
  *
  * @param approx_taylor use Taylor series expansion for the expected value of exp(lambda * c)?
  *
  * @return A vector with the corresponding coefficients of variation
  */
vector cv_dPCR_pre(vector lambda, real nu_pre, vector m, real c, int pre_type, int approx_taylor) {
  if (nu_pre == 0) {
    return cv_dPCR(lambda, m, c);
  }
  int N = num_elements(lambda);
  vector[N] var_exp = nu_pre^2 * lambda^2;
  vector[N] E_exp;
  if (approx_taylor) {
    E_exp = E_exp_taylor(lambda, nu_pre, c);
  } else if (pre_type == 0) {
    E_exp = E_exp_gamma(lambda, nu_pre, c);
  } else if (pre_type == 1) {
    E_exp = E_exp_lnorm(lambda, nu_pre, c);
  } else {
    reject("Unknown pre-PCR noise type");
  }
  vector[N] exp_var = (E_exp - 1) ./ (m * c^2);
  vector[N] total_var = var_exp + exp_var;
  vector[N] cv = sqrt(total_var) ./ lambda;
  return cv;
}

/**
  * Compute the expected value of exp(lambda * c), i.e. evaluate the MGF of the
  * concentration in the PCR at the point t, using a Taylor series expansion
  *
  * Note this only works well for positive t, the series diverges for negative t
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV
  *
  * @param t point of evaluation (moment). Can be a vector.
  *
  * @return A vector with the expected values
  */
vector E_exp_taylor(vector lambda, real nu_pre, vector t) {
  int N = num_elements(lambda);
  vector[N] lambda_t = lambda .* t;
  return(1 + lambda_t + (lambda_t)^2 * (1 + nu_pre^2) / 2);
}
vector E_exp_taylor(vector lambda, real nu_pre, real t) {
  return(E_exp_taylor(lambda, nu_pre, rep_vector(t, num_elements(lambda))));
}

/**
  * Compute the expected value of exp(lambda * c), i.e. evaluate the MGF of the
  * concentration in the PCR at the point c (conversion factor), assuming
  * gamma distributed pre-PCR noise
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV
  *
  * @param t point of evaluation (moment). Can be a vector.
  *
  * @return A vector with the expected values
  */
vector E_exp_gamma(vector lambda, real nu_pre, vector t) {
  real nu_pre2 = nu_pre^2;
  return((1 - lambda .* t * nu_pre2)^(-1/nu_pre2));
}
vector E_exp_gamma(vector lambda, real nu_pre, real t) {
  return(E_exp_gamma(lambda, nu_pre, rep_vector(t, num_elements(lambda))));
}

/**
  * Compute the expected value of exp(lambda * c), i.e. evaluate the MGF of the
  * concentration in the PCR at the point c (conversion factor), assuming
  * gamma distributed pre-PCR noise, on the log scale
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV
  *
  * @param t point of evaluation (moment). Can be a vector.
  *
  * @return A vector with the expected values on the log scale
  */
vector log_E_exp_gamma(vector lambda, real nu_pre, vector t) {
  real nu_pre2 = nu_pre^2;
  return(-log1m(lambda .* t * nu_pre2)/nu_pre2);
}
vector log_E_exp_gamma(vector lambda, real nu_pre, real t) {
  return(log_E_exp_gamma(lambda, nu_pre, rep_vector(t, num_elements(lambda))));
}

/**
  * Compute the expected value of exp(lambda * c), i.e. evaluate the MGF of the
  * concentration in the PCR at the point c (conversion factor), assuming
  * log-normal distributed pre-PCR noise
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV
  *
  * @param t point of evaluation (moment). Can be a vector.
  *
  * @return A vector with the expected values
  */
vector E_exp_lnorm(vector lambda, real nu_pre, vector t) {
  int N = num_elements(lambda);
  vector[N] w = lambert_w0(-t .* lambda * log(1+nu_pre^2) / sqrt(1+nu_pre^2));
  return(exp(-(w^2 + 2*w)/(2*log(1+nu_pre^2))) ./ sqrt(1+w));
}
vector E_exp_lnorm(vector lambda, real nu_pre, real t) {
  return(E_exp_lnorm(lambda, nu_pre, rep_vector(t, num_elements(lambda))));
}

/**
  * Compute the expected value of exp(lambda * c), i.e. evaluate the MGF of the
  * concentration in the PCR at the point c (conversion factor), assuming
  * log-normal distributed pre-PCR noise, on the log scale
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV
  *
  * @param t point of evaluation (moment). Can be a vector.
  *
  * @return A vector with the expected values on the log scale
  */
vector log_E_exp_lnorm(vector lambda, real nu_pre, vector t) {
  int N = num_elements(lambda);
  vector[N] w = lambert_w0(-t .* lambda * log(1+nu_pre^2) / sqrt((1+nu_pre^2)));
  return((-(w^2 + 2*w)/(2*log(1+nu_pre^2))) - log1p(w)/2);
}
vector log_E_exp_lnorm(vector lambda, real nu_pre, real t) {
  return(log_E_exp_lnorm(lambda, nu_pre, rep_vector(t, num_elements(lambda))));
}

// --------------------------------------------------------
// dPCR likelihood based on generalized binomial coefficient
// --------------------------------------------------------

real dPCR_lpdf(real y, real lambda, real m, real c, int norm_l, int norm_u) {
  real p = 1 - exp(-lambda * c);
  real norm_mass_log = 0;
  int norm_n = norm_u - norm_l + 1;
  if (norm_n >= 0) {
    vector[norm_n] all_counts = linspaced_vector(norm_n, norm_l, norm_u);
    vector[norm_n] all_mass = binom_approx_lpdfs(all_counts, m, p);
    norm_mass_log = log_sum_exp(all_mass);
  }
  real obs_mass = binom_approx_lpdf(m * (1 - exp(-y * c)) | m, p);
  return obs_mass - norm_mass_log;
}

real dPCR_lpdf(vector y, vector lambda, vector m, real c, int norm_l, int norm_u) {
  int N = num_elements(y);
  real lp = 0;
  for (i in 1:N) {
    lp += dPCR_lpdf(y[i] | lambda[i], m[i], c, norm_l, norm_u);
  }
  return lp;
}

// vectorized version without normalization
real dPCR_lpdf(vector y, vector lambda, vector m, real c) {
  int N = num_elements(y);
  vector[N] p = 1 - exp(-lambda * c);
  real obs_mass = binom_approx_lpdf(m .* (1 - exp(-y * c)) | m, p);
  return obs_mass;
}

vector dPCR_rng(vector lambda, array[] int m, real c) {
  int N = num_elements(lambda);
  vector[N] p = 1 - exp(-lambda * c);
  array[N] int counts = binomial_rng(m, p);
  vector[N] concs = -log1m(to_vector(counts) ./ to_vector(m)) * (1 / c);
  return concs;
}

// --------------------------------------------------------
// dPCR likelihood based on generalized binomial coefficient
// (conditioned on non-zero measurements)
// --------------------------------------------------------

real dPCR_nonzero_lpdf(vector y, vector lambda, vector m, real c) {
  int N = num_elements(y);
  vector[N] p = 1 - exp(-lambda * c);
  real obs_mass = binom_approx_lpdf(m .* (1 - exp(-y * c)) | m, p);
  // normalize
  obs_mass += -sum(log1m_exp(binom_approx_lpdfs(rep_vector(0.0, N), m, p)));
  return obs_mass;
}

// --------------------------------------------------------
// dPCR likelihood based on generalized binomial coefficient.
// Integration over counts and matching to observed concentrations
// --------------------------------------------------------

/**
  * Squared exponential kernel with magnitude 1.
  *
  * @param y Reference point.
  *
  * @param x Vector with comparison points.
  *
  * @param sigma Length scale parameter.
  *
  * @return A vector with the kernel evaluated for the differences
  * between y and the elements of x.
  */
vector log_se_kernel(real y, vector x, real sigma) {
    return -0.5 * square(y - x) / square(sigma);
  }

/**
  * Normalized squared exponential kernel. This corresponds to a gaussian
  * density conditioned on y>0.
  *
  * @param y Reference point.
  *
  * @param x Vector with comparison points.
  *
  * @param sigma Length scale parameter.
  *
  * @return A vector with the kernel evaluated for the differences
  * between y and the elements of x.
  */
vector log_se_kernel_norm(real y, vector x, real sigma) {
  return -0.5 * square(y - x) / square(sigma)
         - log(sigma)
         - 0.5 * log(2 * pi())
         - log1m(Phi((0 - x) / sigma));
}

/**
  * Integrate over possible integer counts of positive partitions and compare
  * the implied concentration to the observed concentration via a kernel
  *
  * @param y Observed concentration.
  *
  * @param lambda Expected concentration.
  *
  * @param m Number of valid partitions.
  *
  * @param c Conversion factor.
  *
  * @param int_l Lower boundary for number of positive partitions.
  *
  * @param int_u Upper boundary for number of positive partitions.
  *
  * @param sigma Length scale of kernel for comparison between implied
  * concentrations and observed concentrations. Smaller values lead to stricter
  * matching (and slower sampling).
  *
  * @return Probability integrated over all positive partitions counts between
  * int_l and int_u.
  */
real dPCR_int_counts_lpdf(real y, real lambda, real m, real c, int int_l, int int_u, real sigma) {
  real p = 1 - exp(-lambda * c);
  int int_n = int_u - int_l + 1;
  vector[int_n] all_counts = linspaced_vector(int_n, int_l, int_u);
  vector[int_n] all_count_mass = binom_approx_lpdfs(all_counts, m, p);
  vector[int_n] all_match_mass = log_se_kernel_norm(
    y, -1/c * log(soft_lower(1-all_counts/m, 1e-10, 1e10)), sigma
    );
  return log_sum_exp(all_count_mass + all_match_mass);
}

real dPCR_int_counts_lpdf(vector y, vector lambda, vector m, real c, array[] int int_l, array[] int int_u, real sigma) {
  int N = num_elements(y);
  real lp = 0;
  for (i in 1:N) {
    lp += dPCR_int_counts_lpdf(y[i] | lambda[i], m[i], c, int_l[i], int_u[i], sigma);
  }
  // conditioning on non-zero measurements
  lp += -sum(log1m_exp(binom_approx_lpdfs(rep_vector(0.0, N), m, 1 - exp(-lambda * c))));
  return lp;
}
