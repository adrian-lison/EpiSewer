/**
  * Returns a hurdle model probability on the logit scale
  *
  * @param x vector with inputs
  *
  * @param threshold_log the hurdle threshold
  *
  * @param sharpness_log the sharpness of the threshold
  *
  * @return The probability of being below the threshold, on the logit scale
  */
vector hurdle_smooth(vector x, real threshold, real sharpness) {
  return((threshold - x) / threshold * sharpness);
}

/**
  * Returns a hurdle model probability on the logit scale from log scale inputs
  *
  * Uses log_sum_exp-style tricks for better efficiency and avoid direct
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
vector log_hurdle_smooth(vector x_log, real threshold_log, real sharpness_log) {
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
  return(hurdle);
}

