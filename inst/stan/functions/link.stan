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
