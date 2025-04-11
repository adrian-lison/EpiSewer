/**
  * Compute coefficient of variation (CV) in dPCR quantification as a function
  * of the concentration (no pre-PCR noise)
  *
  * @param lambda vector with concentrations (typically gc/mlWW)
  *
  * @param nu_pre pre-PCR CV (all other noise factors)
  *
  * @param m number of partitions
  *
  * @param c conversion factor c = sv, i.e. should be partition volume v
  * multiplied by scaling factor s (concentration in assay / concentration in wastewater)
  *
  * @param n number of replicates that are averaged
  *
  * @return A vector with the corresponding coefficients of variation
  */
vector cv_dPCR(vector lambda, vector m, real c, vector n) {
  int N = num_elements(lambda);
  vector[N] total_var = (exp(lambda * c) - 1) ./ (n .* m * c^2);
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
  * @param m number of partitions
  *
  * @param c conversion factor c = sv, i.e. should be partition volume v
  * multiplied by scaling factor s (concentration in assay / concentration in wastewater)
  *
  * @param n number of replicates that are averaged
  *
  * @param pre_type type of pre-PCR noise (0: gamma, 1: log-normal)
  *
  * @param approx_taylor use Taylor series expansion for the expected value of exp(lambda * c)?
  *
  * @return A vector with the corresponding coefficients of variation
  */
vector cv_dPCR_pre(vector lambda, real nu_pre, vector m, real c, vector n, int pre_type, int approx_taylor) {
  if (nu_pre == 0) {
    return cv_dPCR(lambda, m, c, n);
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
  vector[N] exp_var = (E_exp - 1) ./ (n .* m * c^2);
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

/**
  * Draw number of total valid partitions in dPCR run, assuming
  * lognormally distributed partition loss, and assuming that
  * the maximum number of partitions is 3 standard deviations above the mean
  *
  * @param mu_m Expected number of partitions
  *
  * @param m_cv Coefficient of variation of number of partitions
  *
  * @param noise_raw Standard normal noise for non-centered parameterization
  *
  * @return A vector of total partition numbers, same length as noise_raw
  */
vector total_partitions_noncentered(real m_mu, real m_cv, vector noise_raw) {
  int N = num_elements(noise_raw);
  if (m_cv == 0) {
    return rep_vector(m_mu, N);
  }
  real max_partitions = m_mu + 3 * m_mu * m_cv;
  real lnorm_unit_mean = max_partitions - m_mu;
  real lnorm_unit_sd = m_mu * m_cv;
  real sigma2 = log((lnorm_unit_sd / lnorm_unit_mean)^2 + 1);
  real meanlog = log(lnorm_unit_mean) - sigma2 / 2;
  real sdlog = sqrt(sigma2);
  vector[N] lost_partitions = exp(meanlog + sdlog * noise_raw); // log-normal
  vector[N] total_partitions = max_partitions - lost_partitions;
  total_partitions = softplus(total_partitions, 10); // softplus as soft >0 constraint
  return(total_partitions);
}
