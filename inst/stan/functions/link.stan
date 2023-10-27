/** ---------------------------------------------------------------------------
Helper for link functions
---------------------------------------------------------------------------- */

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
  } else {
    reject("Link function must be one of inv_softplus (0) or scaled_logit (1)");
  }
}
