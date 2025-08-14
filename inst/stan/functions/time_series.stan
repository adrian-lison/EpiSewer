/** ---------------------------------------------------------------------------
Time series process functions
---------------------------------------------------------------------------- */

  /**
  * Vectorized AR(1) process with coefficient of 1
  * Note that if the first value of the process is already included in the
  * increments_standardized vector, then you can provide a start_value = 0
  */
  vector random_walk(vector start_values, vector increments, int diff_order) {
    if (diff_order == 0) {
       return cumulative_sum(append_row(start_values[1], increments));
     } else {
      vector[diff_order] next_start = start_values[2:(diff_order+1)];
      int next_n = num_elements(increments) + diff_order;
      vector[next_n] diffs = random_walk(next_start, increments, diff_order - 1);
      return cumulative_sum(append_row(start_values[1], diffs));
    }
  }

  /**
  * Vectorized MA(q) process with all coefficients = 1
  * Note that if the first value of the process is already included in the
  * increments_standardized vector, then you can provide a start_value = 0
  */
  vector simple_ma(vector start_values, vector increments, int Q, int diff_order) {
    if (diff_order == 0) {
      int n = num_elements(increments)+1;
      vector[n] normal_cumsum = cumulative_sum(append_row(0, increments));
      vector[n] lagged_cumsum = cumulative_sum(append_row(rep_vector(0, Q+2), increments))[1:n];
      vector[n] window_sum = normal_cumsum - lagged_cumsum;
      return start_values[1] + window_sum; // start_value represents mean here
     } else {
      vector[diff_order] next_start = start_values[2:(diff_order+1)];
      int next_n = num_elements(increments) + diff_order - 1;
      vector[next_n] diffs = simple_ma(next_start, increments, Q, diff_order - 1);
      return cumulative_sum(append_row(start_values[1], diffs));
    }
  }

  /**
  * Damped Holt's method as an innovation state space model, vectorized implementation
  * Linear method with additive noise and without seasonality
  */
  vector holt_damped_process(vector start_values, real alpha, real beta_star, real phi, vector noise, int diff_order) {
    if (diff_order == 0) {
      int n = 1 + num_elements(noise); // n = 1 (for start value) + length of noise
      vector[n] y; // observations (from 1:n)
      vector[n] epsilons = append_row(0, noise); // innovations (from 1:n)
      // We want y_1 to be equal to the start value/intercept, hence epsilon_1 must have zero noise.
      // --> epsilons = [epsilon_1, epsilon_2, ...] = [0, noise_1, ...]
      // - 2) We also need to model epsilon_0, as we need it for l_0 and b_0 (they define y_1). Since we neither want epsilon_0 to affect y_1, we also set it to zero.
      // Therefore, we add two zeros to the head of the noise vector:


      // Level values (from l_0:l_n-1)
      // --> l_{t} = l_0 + sum(epsilon_{1:t})
      vector[n] level;
      level[1] = start_values[1]; // l_0
      level[2:n] = start_values[1] + alpha * cumulative_sum(epsilons[1:(n-1)]); // l_1:l_n-1, note index shift for epsilon

      // Trend values (from b_0:b_n-1)
      vector[n] trend;
      real beta = alpha * beta_star; // below, we always need the product of alpha and beta*
      if (phi == 0) {
        // special case: no trend
        trend = rep_vector(0, n);
      } else {
        trend[1] = 0; // b_0, assume no trend for y_0 -> y_1 (so that y_1 corresponds to the supplied intercept without trend)
        trend[2] = start_values[2]; // b_1, assume that trend starts for step y_1 -> y_2
        if (phi == 1) {
          // special case: trend, no dampening (vectorized)
          trend[3:n] = trend[2] + beta * cumulative_sum(epsilons[2:(n-1)]); // b_2:b_n-1, note index shift for epsilon
        } else {
          // general case: trend and dampening
          for (t in 3:n) { trend[t] = phi * trend[t - 1] + beta * epsilons[t-1]; } // b_2:b_n-1, note index shift for epsilon
        }
        // update level: l_t = l_{t-1} + b{t-1} + alpha * epsilon_t
        // --> l_t = alpha * sum(epsilon_{0:t}) + sum(b_{0:t-1})
        // --> For l_0, we also need b_{-1} (assume it is 0)!
        level = level + phi * cumulative_sum(append_row(0,trend)[1:n]); // update l_0:l_n-1
      }

      // Compute observations
      // the y values correspond to y_t, i.e. from y_1 to y_T
      // y_t = l_{t-1} + b_{t-1} + epsilon_t
      y = level + phi * trend + epsilons[1:n];
      return y;
    } else {
      vector[diff_order+1] next_start = start_values[2:(diff_order+2)];
      int next_n = num_elements(noise) + diff_order;
      // Compute holt_damped_process at one diff_order lower
      vector[next_n] diffs = holt_damped_process(next_start, alpha, beta_star, phi, noise, diff_order - 1);
      // Add up differences
      return cumulative_sum(append_row(start_values[1], diffs));
    }
  }

  vector ar1_process(real phi, vector noise) {
    int n = num_elements(noise);
    vector[n] y; // observations (from 1:n)
    y[1] = noise[1];
    for (t in 2:n) {
      y[t] = phi * y[t-1] + noise[t];
    }
    return(y);
  }


  vector ar1_process(real start_value, real phi, vector noise) {
    int n = num_elements(noise) + 1; // n = 1 (for start value) + length of noise
    vector[n] y; // observations (from 1:n)
    vector[n] epsilons = append_row(0, noise); // innovations (from 1:n)
    y[1] = start_value;
    for (t in 2:n) {
      y[t] = phi * y[t-1] + epsilons[t];
    }
    return(y);
  }

  real ets_coefficient_priors_lp(
    array[] real ets_alpha, array[] real ets_alpha_prior,
    array[] real ets_beta, array[] real ets_beta_prior,
    array[] real ets_phi, array[] real ets_phi_prior
  ) {
    real tar = 0;
    tar += beta2_prior_lpdf(ets_alpha | ets_alpha_prior);
    tar += beta2_prior_lpdf(ets_beta | ets_beta_prior);
    tar += beta2_prior_lpdf(ets_phi | ets_phi_prior); // dampening needs a tight prior, roughly between 0.8 and 0.98
    return(tar);
  }

  /**
  * Dampen the trend of a time series
  *
  * @param start The start value of the time series. This is used to compute
  *   the trend of the first value.
  * @param values The actual values of the time series that should be dampened.
  * @param dampening The exponential dampening parameter to apply. A value
  *   of 1 means no dampening, a value of 0 means no trend (flat time series).
  *
  * @return A vector of the dampened time series, same length as `values`.
  */
  vector dampen_trend(real start, vector values, real dampening) {
    int n = num_elements(values);
    vector[n] dampened_trend;
    dampened_trend[1] = dampening * (values[1] - start);
    for (i in 2:n) {
      dampened_trend[i] = dampening^i * (values[i]-values[i-1]);
    }
    return(start + cumulative_sum(dampened_trend));
  }
