real gev_lpdf(vector y, real mu, real sigma, real xi) {
    int N = num_elements(y);
    vector[N] z = 1 + (y - mu) * xi / sigma;
    vector[N] lp = (1 + (1 / xi)) * log(z) + pow(z, -1/xi);
    return -sum(lp) - N * log(sigma);
}

vector gev_rng(vector mu, real sigma, real xi) {
    int N = num_elements(mu);
    vector[N] u = to_vector(uniform_rng(rep_vector(0, N), 1));
    if (xi != 0) {
      return mu + sigma * (pow(-log(u), -xi) - 1) / xi;
    } else {
      return mu - sigma * log(-log(u));
    }
}

