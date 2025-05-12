/** ---------------------------------------------------------------------------
Helper for link functions
---------------------------------------------------------------------------- */

  /**
  * Softplus activation function
  */
  real softplus(real x, real k) {
    return(log1p_exp(k * x) / k);
  }

  /**
  * Softplus activation function
  */
  vector softplus(vector x, real k) {
    return(log1p_exp(k * x) / k);
  }

  /**
  * Soft upper limit using the softplus function
  */
  real soft_upper(real x, real u, real k) {
    return(u - softplus(u - x, k));
  }

  /**
  * Soft upper limit using the softplus function
  */
  vector soft_upper(vector x, real u, real k) {
    return(u - softplus(u - x, k));
  }

  /**
  * Soft lower limit using the softplus function
  */
  real soft_lower(real x, real l, real k) {
    return(l + softplus(x - l, k));
  }

  /**
  * Soft lower limit using the softplus function
  */
  vector soft_lower(vector x, real l, real k) {
    return(l + softplus(x - l, k));
  }

  /**
  * Logistic
  */
  real logistic(real x, real c, real a, real k) {
    return(c / (1 + a * exp(-k * x)));
  }

  /**
  * Logistic
  */
  vector logistic(vector x, real c, real a, real k) {
    return(c / (1 + a * exp(-k * x)));
  }

    /**
  * Helper function to compute the softmax-weighted average of a vector, also
  * known as Boltzman operator.
  *
  * @param x Input vector.
  *
  * @param beta Temperature parameter. The operator approaches the maximum operator
  *  when beta approaches infinity, the mean when beta approaches zero, and the
  *  minimum operator when beta approaches minus infinity.
  *
  * @return The softmax-weighted average.
  */
  real boltzmann(vector x, real beta) {
    int n = num_elements(x);
    vector[n] weights = exp(beta * x);
    return dot_product(weights / sum(weights), x);
  }

  vector cumulative_boltzmann(vector x, real beta) {
    int n = num_elements(x);
    vector[n] cum_boltzmann = rep_vector(0, n);
    vector[n] weights = exp(beta * x);
    for (i in 1:n) {
      cum_boltzmann[i] = dot_product(weights[1:i]/sum(weights[1:i]), x[1:i]);
    }
    return(cum_boltzmann);
  }

  vector boltzmann_elementwise(vector x, vector y, real beta) {
    if (num_elements(x) != num_elements(y)) {
      reject("x and y must have the same number of elements");
    }
    int n = num_elements(x);
    vector[n] cum_boltzmann;
    vector[n] weights_x = exp(beta * x);
    vector[n] weights_y = exp(beta * y);
    vector[n] weights_sum = weights_x + weights_y;
    return(weights_x ./ weights_sum .* x + weights_y ./ weights_sum .* y);
  }

  /**
  * Apply link function to a vector
  */
  vector apply_link(vector x, array[] real link) {
    if (link[1] == 0) {
      return(softplus(
        x,
        link[2] // hyperparameter k (sharpness)
        ));
    } else if (link[1] == 1) {
      return(logistic(
        x,
        link[2], // hyperparameter c (maximum R)
        link[3], // hyperparameter a
        link[4] // hyperparameter k
        ));
    } else if (link[1] == 2) {
      return(exp(x));
    }
    else {
      reject("Link function must be one of inv_softplus (0), scaled_logit (1), or log (2)");
    }
  }
