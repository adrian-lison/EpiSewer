/**
  * Returns log probability of a logit sigmoid hurdle model
  *
  * The model uses a sigmoid function on the logit scale which has its
  * inflection point at a specified threshold.
  *
  * @param x vector with inputs
  *
  * @param threshold_log the hurdle threshold
  *
  * @param sharpness_log the sharpness of the threshold
  *
  * @return The probability of being below the threshold, on the logit scale
  */
vector log_hurdle_sigmoid(vector x, real threshold, real sharpness) {
  return(log_inv_logit((threshold - x) / threshold * sharpness));
}

/**
  * Returns log probability of a logit sigmoid hurdle model
  * from log scale parameters
  *
  * The model uses a sigmoid function on the logit scale which has its
  * inflection point at a specified threshold.
  *
  * Uses log_sum_exp-style tricks for better efficiency and avoids direct
  * exponentiation of the inputs. This is numerically more stable, especially
  * when the input and the threshold are not too far apart.
  *
  * @param x vector with inputs (log scale)
  *
  * @param threshold_log the hurdle threshold (log scale)
  *
  * @param sharpness_log the sharpness of the threshold (log scale)
  *
  * @return The probability of being below the threshold, on the logit scale
  */
vector log_hurdle_sigmoid_log(vector x_log, real threshold_log, real sharpness_log) {
  int n = num_elements(x_log);
  vector[n] hurdle = rep_vector(0, n);
  for (i in 1:n) {
    if (threshold_log > x_log[i]) {
      hurdle[i] = exp(
        log1m_exp(x_log[i] - threshold_log) + sharpness_log
      );
    } else if (threshold_log < x_log[i]) {
      hurdle[i] = -exp(
        x_log[i] + log1m_exp(threshold_log - x_log[i]) - threshold_log + sharpness_log
      );
    }
  }
  return(log_inv_logit(hurdle));
}

/**
  * Returns log probability of an exponential hurdle model assuming Gamma-distributed noise
  *
  * The model uses an exponential function, the probability becomes 1 as x goes to zero.
  * Moreover, gamma-distributed noise around the inputs can be modeled. The
  * stronger the noise is, the longer the probability will stay close to 1 even
  * for larger values of x.
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
  * @return The probability of being below the threshold, on the logit scale
  */
vector log_hurdle_exponential(vector lambda, vector hurdle_scale, real nu_pre, int pre_type) {
  if (nu_pre == 0) {
    return(lambda .* (-hurdle_scale));
  } else {
    if (pre_type == 0) {
        return(log_E_exp_gamma(lambda, nu_pre, -hurdle_scale));
    } else if (pre_type == 1) {
        return(log_E_exp_lnorm(lambda, nu_pre, -hurdle_scale));
    } else {
    reject("Unknown pre-PCR noise type");
    }
  }
}
