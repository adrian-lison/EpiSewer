  real approx_poisson_log_lpdf(vector y, vector mean_log) {
    int n = num_elements(y);
    vector[n] sigma2 = log1p_exp(-mean_log);
    vector[n] mu = mean_log - sigma2/2;
    return normal_lpdf(y | mu, sqrt(sigma2));
  }
  
  real approx_negative_binomial_log_lpdf(vector y, vector mean_log, real xi) {
    int n = num_elements(y);
    //variance = mean * (1 + mean*square(xi)) = mean + square(mean)*square(xi)
    //sigma2 = log(1+variance/square(mean)) = log(1 + 1/mean + square(xi))
    vector[n] sigma2 = log1p_exp(log_sum_exp_elementwise(-mean_log, rep_vector(2*log(xi), n)));
    vector[n] mu = mean_log - sigma2/2;
    return normal_lpdf(y | mu, sqrt(sigma2));
  }
  
  real approx_negative_binomial_log_noncentered(real mean_log, real xi, real raw_y) {
    real sigma2 = log1p_exp(log_sum_exp(-mean_log, 2*log(xi)));
    real mu = mean_log - sigma2/2;
    return mu + raw_y * sqrt(sigma2);
  }
  
  vector approx_negative_binomial_log_noncentered(vector mean_log, real xi, vector raw_y) {
    int n = num_elements(raw_y);
    vector[n] sigma2 = log1p_exp(log_sum_exp_elementwise(-mean_log, rep_vector(2*log(xi), n)));
    vector[n] mu = mean_log - sigma2/2;
    return mu + raw_y .* sqrt(sigma2);
  }
