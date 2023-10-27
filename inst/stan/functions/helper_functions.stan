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
  * Softplus activation function
  */
  vector softplus(vector x, real k) {
    return(log1p_exp(k * x) / k);
  }

  /**
  * Logistic
  */
  vector logistic(vector x, real c, real a, real k) {
    return(c / (1 + a * exp(-k * x)));
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
