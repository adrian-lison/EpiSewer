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
* @param break_dist Default distance between breakpoints. The segments all have
*   the same length by default (equidistant breakpoints).
*
* @param min_break_dist Minimum distance between breakpoints. This will
*   implement a constraint on how close to breakpoints can be to each other.
*
* @param delays Array of vectors. Each vector has length break_dist, and there
*   is one for each segment. The vector elements indicate how much the default
*   breakpoint position gets delayed, where 0 means "place breakpoint here" and
*   1 means "delay breakpoint further"
*
* @param logistic_k Sharpness parameter of logistic function. This influences
*   how strict the approximation to the piecewise constant function is. Higher
*   values mean sharper jumps / step-like behaviour but can slow down sampling.
*
* @param skip_tol Tolerance for `min_break_dist` constraint. When two segments
*  differ by less than `skip_tol`, they are considered identical, meaning that
*. the next segment can is not constrained by `min_break_dist` anymore.
*/
vector soft_changepoint(real intercept, vector segments,
                        int length_intercept, int length_total,
                        int break_dist, int min_break_dist,
                        array[] vector delays, real logistic_k,
                        real skip_tol, real skip_tol_logistic_k) {
  int n_segments = num_elements(segments);
  array[n_segments] vector[break_dist] break_delay;

  // knot placement
  for (seg_i in 1:n_segments) {

    // check if segment qualifies as skip
    real skip = is_skip(
      (seg_i > 1 ? segments[seg_i-1] : intercept), segments[seg_i],
      skip_tol, skip_tol_logistic_k
      );

    // compute breakpoint delays
    break_delay[seg_i] = 1-cumulative_sum(delays[seg_i]);

    // apply min_break_dist constraint
    if ((seg_i > 1) && (min_break_dist > 1)) {
      break_delay[seg_i][1:(min_break_dist-1)] += (1-skip) * break_delay[seg_i - 1][
        (break_dist-min_break_dist+1):(break_dist-1)
        ];
    }

    // logistic link: ensure that break_delay is in [0,1], and increase sharpness
    break_delay[seg_i] = 1 / (1 + exp(logistic_k * (0.5 - break_delay[seg_i])));
  }

  // compute function values
  vector[length_total] y;
  if (length_intercept>=1) {
    y[1:length_intercept] = rep_vector(intercept, length_intercept);
  }
  y[(length_intercept+1):(length_intercept+break_dist)] = (
    intercept * break_delay[1] + segments[1] * (1-break_delay[1])
    );
  for (seg_i in 2:n_segments) {
    int seg_start = length_intercept + 1 + (seg_i-1) * break_dist;
    int seg_end = length_intercept + seg_i * break_dist;
    y[seg_start:seg_end] = (
      segments[seg_i - 1] * break_delay[seg_i] +
      segments[seg_i] * (1-break_delay[seg_i])
      );
  }
  int last_knot = length_intercept+(n_segments*break_dist);
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
