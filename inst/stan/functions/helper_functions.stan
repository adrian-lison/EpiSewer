/** ---------------------------------------------------------------------------
Helper functions for primitive operations
---------------------------------------------------------------------------- */

  /**
  * Cumulative product for column vector
  */
  vector cumulative_prod_column(vector x) {
    return exp(cumulative_sum(log(x)));
  }

  /**
  * Cumulative product of 1-x for column vector
  */
  vector cumulative_prod1m_column(vector x) {
    return exp(cumulative_sum(log1m(x)));
  }

  /**
  * Cumulative product for row vector
  */
  row_vector cumulative_prod_row(row_vector x) {
    return exp(cumulative_sum(log(x)));
  }

  /**
  * Cumulative product of 1-x for row vector
  */
  row_vector cumulative_prod1m_row(row_vector x) {
    return exp(cumulative_sum(log1m(x)));
  }

  /**
  * Element-wise log_sum_exp between two vectors
  */
  vector log_sum_exp_elementwise(vector x, vector y) {
    int n = num_elements(x); // we assume x and y have equal length
    vector[n] res;
    for (i in 1:n) {
      res[i] = log_sum_exp(x[i], y[i]);
    }
    return(res);
  }

  /**
  * Element-wise log_diff_exp between two vectors
  */
  vector log_diff_exp_elementwise(vector x, vector y) {
    int n = num_elements(x); // we assume x and y have equal length
    vector[n] res;
    for (i in 1:n) {
      res[i] = log_diff_exp(x[i], y[i]);
    }
    return(res);
  }

  /**
  * Efficient dot product on log scale
  */
  real log_dot_product(vector x, vector y) {
    return(log_sum_exp(x + y));
  }

  /**
  * Efficient cumulative sum on log scale
  */
  vector log_cumulative_sum(vector x) {
    int n = num_elements(x);
    vector[n] cx;
    cx[1] = x[1];
    for (i in 2:n) {
      cx[i] = log_sum_exp(cx[i-1], x[i]);
    }
    return(cx);
  }

  /**
  * Vectorized rolling sum, right-aligned
  */
  vector rolling_sum_r(vector x, int window_size) {
    if (window_size == 1) {
      return(x);
    } else {
      int n = num_elements(x);
      vector[n] cx = cumulative_sum(x);
      vector[n] cx_shifted;
      // cx_shifted[1:(window_size-1)] left as NaN
      cx_shifted[window_size] = 0;
      cx_shifted[(window_size+1):n] = cx[1:(n-window_size)];
      return(cx - cx_shifted);
    }
  }

  /**
  * Vectorized rolling sum, right-aligned, on log scale
  */
  vector log_rolling_sum_r(vector x, int window_size) {
    if (window_size == 1) {
      return(x);
    } else {
      int n = num_elements(x);
      vector[n] cx = log_cumulative_sum(x);
      vector[n] cx_shifted;
      // cx_shifted[1:(window_size-1)] left as NaN
      cx_shifted[window_size] = 0;
      cx_shifted[(window_size+1):n] = cx[1:(n-window_size)];
      return(log_diff_exp_elementwise(cx, cx_shifted));
    }
  }

  /**
  * Vectorized rolling mean, right-aligned
  */
  vector rolling_mean_r(vector x, int window_size) {
    return(rolling_sum_r(x, window_size) / window_size);
  }

  /**
  * Vectorized rolling mean, right-aligned, on log scale
  */
  vector log_rolling_mean_r(vector x, int window_size) {
    return(log_rolling_sum_r(x, window_size) - log(window_size));
  }

  /**
  * Count number of zero entries in an array
  */
  int num_zeros(array[] int y) {
    int sum = 0;
    for (n in 1:size(y)) {
      sum += (y[n] == 0);
    }
    return sum;
  }

  /**
  * Count number of zero entries in a vector
  */
  int num_zeros(vector y) {
    int sum = 0;
    for (n in 1:num_elements(y)) {
      sum += (y[n] == 0);
    }
    return sum;
  }

  /**
  * Return a parameter - or a fixed value if its prior has zero variance
  *
  * @param param The parameter. This assumes that the parameter is stored inside
  * an array of length 1 if the parameter is not fixed and length 0 otherwise
  *
  * @param prior The prior for the parameter. This assumes that the prior is
  * stored in an array of length 2, where the first element contains the mean
  * and the second element the standard deviation of the prior.
  *
  * @return If the prior has nonzero variance, the parameter is returned, if
  * instead the prior has zero variance, then the mean of the prior is returned.
  */
  real param_or_fixed(array[] real param, array[] real prior) {
    if (prior[2] > 0) {
      return(param[1]);
    } else {
      return(prior[1]);
    }
  }

  /**
  * Check if a vector contains NaN values
  */
  int has_nan_vector(vector x) {
    int has_nan = 0;
    for (i in 1:num_elements(x)) {
      if (is_nan(x[i])) {
        has_nan = 1;
      }
    }
    return has_nan;
  }

  /*
  Function to get the indices of NaN values in a vector
  */
  array[] int get_nan_positions(vector x) {
    int nan_count = 0;
    array[num_elements(x)] int positions;

    for (i in 1:num_elements(x)) {
      if (is_nan(x[i])) {
        nan_count += 1;
        positions[nan_count] = i;
      }
    }
    return positions[1:nan_count];
  }

  vector trim_or_reject_lb(vector x, real lb_trim, real lb_reject) {
    int n = num_elements(x);
    array[n] int rej_positions;
    int rej_count = 0;
    for (i in 1:n) {
      if (x[i] < lb_reject) {
        rej_count += 1;
        rej_positions[rej_count] = i;
      }
    }
    if (rej_count > 0) {
      reject(
        "The following vector elements were below ",
        "the lower bound for rejection (", lb_reject, "). ",
        "Indices: ", rej_positions[1:rej_count], " | ",
        "Values: ", x[rej_positions[1:rej_count]]
        );
    } else {
      return(fmax(x, rep_vector(lb_trim, n)));
    }
  }

  vector trim_or_reject_ub(vector x, real ub_trim, real ub_reject) {
    int n = num_elements(x);
    array[n] int rej_positions;
    int rej_count = 0;
    for (i in 1:n) {
      if (x[i] > ub_reject) {
        rej_count += 1;
        rej_positions[rej_count] = i;
      }
    }
    if (rej_count > 0) {
      reject(
        "The following vector elements were above ",
        "the upper bound for rejection (", ub_reject, "). ",
        "Indices: ", rej_positions[1:rej_count], " | ",
        "Values: ", x[rej_positions[1:rej_count]]
        );
    } else {
      return(fmin(x, rep_vector(ub_trim, n)));
    }
  }

  /*
  append_row extended to three elements
  */
  vector append_row3(vector x, vector y, vector z) {
    return append_row(append_row(x, y), z);
  }
