holt_damped_process_noncentered <- function(
    alpha, beta_star, phi, l_start, b_start, increments) {
  increments <- c(0, increments)
  n <- length(increments)
  beta <- alpha * beta_star

  if (phi == 0) {
    sum_b <- rep(0, n)
  } else if (phi == 1) {
    b <- b_start + beta * cumsum(c(0, increments))[1:n]
    sum_b <- cumsum(b)
  } else {
    b <- rep(NA, n)
    b[1] <- b_start
    for (t in 2:n) {
      b[t] <- phi * b[t - 1] + beta * increments[t - 1]
    }
    sum_b <- phi * cumsum(b)
  }

  y <- l_start + alpha * cumsum(c(0, increments))[1:n] + sum_b + increments

  return(y)
}
