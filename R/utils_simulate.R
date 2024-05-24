holt_damped_process <- function(start_values, alpha, beta_star, phi, noise, diff_order) {
  if (diff_order == 0) {
    n <- 1 + length(noise)
    y <- numeric(n)
    epsilons <- c(0, noise)

    level <- numeric(n)
    level[1] <- start_values[1]
    level[2:n] <- start_values[1] + alpha * cumsum(epsilons[1:(n-1)])

    trend <- numeric(n)
    beta <- alpha * beta_star
    if (phi == 0) {
      trend <- rep(0, n)
    } else {
      trend[1] <- 0
      trend[2] <- start_values[2]
      if (phi == 1) {
        trend[3:n] <- trend[2] + beta * cumsum(epsilons[2:(n-1)])
      } else {
        for (t in 3:n) {
          trend[t] <- phi * trend[t - 1] + beta * epsilons[t - 1]
        }
      }
      level <- level + phi * cumsum(c(0, trend)[1:n])
    }

    y <- level + phi * trend + epsilons[1:n]
    return(y)
  } else {
    next_start <- start_values[2:(diff_order + 2)]
    next_n <- length(noise) + diff_order
    diffs <- holt_damped_process(next_start, alpha, beta_star, phi, noise, diff_order - 1)
    return(cumsum(c(start_values[1], diffs)))
  }
}
