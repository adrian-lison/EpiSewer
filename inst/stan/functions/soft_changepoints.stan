/**
* Soft changepoint model approximating a piecewise constant function
*
* @param intercept Intercept / start value of the function.
*
* @param segments Values for the segments of the function.
*
* @param length_intercept Point of the first breakpoint.
*
* @param length_total Total length, e.g. number of time points in a time series.
*
* @param distance Default distance between breakpoints. The segments all have
*   the same length by default (equidistant breakpoints).
*
* @param min_distance Minimum distance between breakpoints. This will
*   implement a constraint on how close to breakpoints can be to each other.
*
* @param delays Array of vectors. Each vector has length distance, and there
*   is one for each segment. The vector elements indicate how much the default
*   breakpoint position gets delayed, where 0 means "place breakpoint here" and
*   1 means "delay breakpoint further"
*
* @param boltzmann_sharpness Sharpness parameter (temperature) for the
* Boltzmann operator (softmax-weighted average). Influences how sharp the
* approximation of the maximum function for the fuzzy OR is.
*
* @param skip_tol Tolerance for `min_distance` constraint. When two segments
*  differ by less than `skip_tol`, they are considered identical, meaning that
*. the next segment can is not constrained by `min_distance` anymore.
*/
vector soft_changepoint(real intercept, vector segments,
                        int length_intercept, int length_total,
                        int distance, int min_distance,
                        array[] vector delays,
                        real boltzmann_sharpness,
                        real skip_tol, real skip_tol_logistic_k) {
  int n_segments = num_elements(segments);
  array[n_segments] vector[distance] break_delay;

  // knot placement
  for (seg_i in 1:n_segments) {

    // check if segment qualifies as skip
    real skip = is_skip(
      (seg_i > 1 ? segments[seg_i-1] : intercept), segments[seg_i],
      skip_tol, skip_tol_logistic_k
      );

    // compute breakpoint delays
    break_delay[seg_i] = 1-cumulative_sum(delays[seg_i]);

    // apply min_distance constraint
    if ((seg_i > 1) && (min_distance > 1)) {
      break_delay[seg_i][1:(min_distance-1)] = boltzmann_elementwise(
        break_delay[seg_i][1:(min_distance-1)],
        (1-skip) * break_delay[seg_i - 1][(distance-min_distance+1):(distance-1)],
        boltzmann_sharpness
      );
    }
  }

  // compute function values
  vector[length_total] y;
  if (length_intercept>=1) {
    y[1:length_intercept] = rep_vector(intercept, length_intercept);
  }
  y[(length_intercept+1):(length_intercept+distance)] = (
    intercept * break_delay[1] + segments[1] * (1-break_delay[1])
    );
  for (seg_i in 2:n_segments) {
    int seg_start = length_intercept + 1 + (seg_i-1) * distance;
    int seg_end = length_intercept + seg_i * distance;
    y[seg_start:seg_end] = (
      segments[seg_i - 1] * break_delay[seg_i] +
      segments[seg_i] * (1-break_delay[seg_i])
      );
  }
  int last_knot = length_intercept+(n_segments*distance);
  y[(last_knot+1):length_total] = rep_vector(
    y[last_knot], length_total-last_knot
    );
  return(y);
}

/**
* Helper function to check if two segments are approximately identical, i.e.
*  the breakpoint has been skipped.
*
* @param x Value of first segment.
*
* @param y Value of second segment.
*
* @param skip_tol Tolerance for being considered identical. When two segments
*  differ by less than `skip_tol`, they are considered identical.
*
* @param logistic_k Sharpness parameter of logistic function to implement
*  `skip_tol` as a soft constraint. The higher `logistic_k`, the stricter the
*. constraint.
*
* @return This function gives a real-valued return that approximates a boolean
*  result. If the segments differ by less than skip_tol, a value close to 1 will
*  be returned, if they differ by more than skip_tol, a value close to 0 will be
*  returned. If the difference is close to skip_tol, however, values in between
*  may be returned as well.
*/
real is_skip(real x, real y, real skip_tol, real logistic_k) {
  if (skip_tol > 0) {
    return(1 - 1 / (1 + exp(logistic_k * (1 - (((x-y)/skip_tol)^2)))));
  } else {
    return(0);
  }
}

/**
* Compute adjusted prior for changepoint model under min_distance constraint
*
* @description By adjusting the prior for the changepoint position, we ensure
* that changepoints occur with equal prior probability at each time step even
* after accounting for the min_distance constraint.
*
* @param d Default distance between breakpoints. The segments all have
*   the same length by default (equidistant breakpoints).
*
* @param d_min Minimum distance between breakpoints. This will
*   implement a constraint on how close to breakpoints can be to each other.
*
* @return A simplex with adjusted weights for the changepoint position prior.
*/
vector changepoint_prior_adjusted(int d, real d_min) {
  vector[d] weights = rep_vector(0, d);
  for (j in 1:d) {
    if (j <= d_min) {
      weights[j] = (d - d_min) / ((d-d_min+j-1) * (d-d_min+j));
    } else {
      weights[j] = (1 - d_min/d) / (d-d_min);
    }
  }
  return(weights);
}
