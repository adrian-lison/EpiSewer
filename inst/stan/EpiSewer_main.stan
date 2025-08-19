functions {
  #include functions/helper_functions.stan
  #include functions/link.stan
  #include functions/time_series.stan
  #include functions/approx_count_dist.stan
  #include functions/renewal.stan
  #include functions/dist_normal.stan
  #include functions/dist_lognormal.stan
  #include functions/dist_gamma.stan
  #include functions/dist_loggamma.stan
  #include functions/dist_gamma_sum.stan
  #include functions/dist_beta.stan
  #include functions/dist_gev.stan
  #include functions/dist_dirichlet.stan
  #include functions/discretize.stan
  #include functions/hurdle.stan
  #include functions/pcr_noise.stan
  #include functions/soft_changepoints.stan
  #include functions/gaussian_process.stan
}
data {
  // Measurements ----
  int<lower=0> T; // number of total days in measured time span
  int<lower=0> h; // number of days to forecast
  int<lower=0> n_samples; // number of samples from different dates
  int<lower=0> n_measured; // number of different measurements (including replicates)
  array[n_samples] int<lower=1, upper=T> sample_to_date; // index mapping samples to dates
  array[n_measured] int<lower=1, upper=n_samples> measure_to_sample; // index mapping measurements to samples
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  vector<lower=0>[n_measured] n_averaged; // number of averaged technical replicates per measurement (is vector for vectorization)
  vector<lower=0>[n_measured] dPCR_total_partitions; // total number of partitions in dPCR
  int<lower=1> w; // composite window: how many past days the samples cover, e.g. 1 for individual day samples, 7 for weekly composite samples, ...
  int<lower=0, upper=3> obs_dist; // Parametric distribution for observation likelihood: 0 (default) for gamma, 1 for log-normal, 2 for truncated normal, 3 for normal

  // Flow ----
  vector<lower=0>[T+h] flow; // flow rate for normalization of measurements

  // Sample date effects model ----
  int<lower=0> K; // number of sample date predictors
  matrix[T+h, K] X; // sample date predictor design matrix
  array[K > 0 ? 2 : 0] real eta_prior; // prior for sample date effects

  // Pre-replicate noise ----
  int<lower=0, upper=1> pr_noise; // Pre-replicate noise: Model variation before replication stage?
  array[pr_noise ? 2 : 0] real nu_psi_prior; // Prior on coefficient of variation for pre-replicate noise

  // Coefficient of variation (CV) of measurements ----
  int<lower=0, upper=2> cv_type; // 0 for constant, 1 for dPCR, 2 for constant_var
  array[2] real nu_upsilon_a_prior; // prior for pre-PCR CV
  int<lower=0, upper=1> total_partitions_observe; // 0 for not observed, 1 for observed
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_b_mu_prior; // prior for parameter 2 of CV formula (number of partitions). Scaled by 1e-4 for numerical efficiency.
  array[cv_type == 1 && total_partitions_observe!=1 ? 2 : 0] real nu_upsilon_b_cv_prior; // prior on partition number variation (coefficient of variation)
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_c_prior; // prior for parameter 3 of CV formula (partition size*(scaling factor, i.e. exp_conc_assay/exp_conc_ww)). Scaled by 1e+5 for numerical efficiency.
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_type; // 0 for gamma, 1 for log-normal
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_approx_taylor; // 0 for no Taylor expansion approximation, 1 for Taylor expansion approximation

  // Limit of detection ----
  // LOD_model = 0: no LOD
  // LOD_model = 1: assumed LOD, LOD_scale provided
  // LOD_model = 2: estimated LOD based on dPCR model, needs dPCR parameters
  int<lower=0, upper=2> LOD_model;
  array[(LOD_model == 1) ? 1 : 0] real<lower=0> LOD_scale;
  real<lower=0, upper=1> LOD_drop_prob; // probability threshold for non-detection below which log likelihood contributions of observed concentrations are dropped from LOD model

  // Residence time ----
  int<lower=0> D; // maximum residence time
  vector[D + 1] residence_dist; // residence time distribution
  // --> probability for residence of zero days (same day arrival at sampling site) comes first

  // Shedding ----
  int<lower=1> S; // maximum number of days with shedding
  int<lower=0, upper=3> shedding_dist_type; // 0 for fixed, non-parametric, 1 for exponential, 2 for gamma, 3 for log-normal
  int<lower=1> shedding_dist_n; // number of shedding load priors
  vector[S + 1] shedding_dist; // Non-parametric shedding load distribution. If shedding_dist_type is parametric and estimated, this is the discretized mean of the prior.
  // --> probability for shedding today comes first
  vector[(shedding_dist_type > 0 && shedding_dist_n > 1) ? shedding_dist_n : 0] shedding_dist_weights_prior; // weights prior for shedding distribution priors
  array[shedding_dist_type > 0 ? shedding_dist_n : 0, shedding_dist_type > 0 ? 2 : 0] real shedding_dist_mean_prior; // priors for mean of shedding distribution
  array[shedding_dist_type > 1 ? shedding_dist_n : 0, shedding_dist_type > 1 ? 2 : 0] real shedding_dist_cv_prior; // priors for cv of shedding distribution

  real<lower=0> load_mean; // mean load shed per person
  int<lower=0, upper=1> load_vari; // model individual-level variation in shedding loads?
  array[load_vari ? 2 : 0] real nu_zeta_prior; // prior on coefficient of variation of individual-level load
  int<lower = 0> n_zeta_normal_approx;
  int<lower = 0> n_zeta_exact;
  array[load_vari ? n_zeta_normal_approx : 0] int<lower = 1, upper = S + D + T> zeta_normal_approx; // dates on which zeta can be approximated
  array[load_vari ? n_zeta_exact : 0] int<lower = 1, upper = S + D + T> zeta_exact; // dates on which zeta should not be approximated

  // Outliers ----
  int<lower=0, upper=1> outliers;
  array[outliers ? 3 : 0] real<lower=0> epsilon_prior;

  // Incubation period ----
  int<lower=1> L; // maximum incubation period
  vector[L + 1] incubation_dist; // incubation period distribution
  // --> probability for a delay of zero comes first

  // Generation interval ----
  int<lower=1> G; // maximum generation interval
  vector<lower=0>[G] generation_dist; // generation interval distribution
  // --> probability for a delay of one comes first (zero excluded)

  // Seeding of infections ----
  int<lower=0, upper=2> seeding_model; // 0 for fixed, 1 for random walk, 2 for growth rate seeding
  array[2] real iota_log_seed_intercept_prior;
  array[seeding_model > 0 ? 2 : 0] real iota_log_seed_sd_prior;
  int<lower=0> se; // seeding extension (used when time series starts with many non-detects)

  // Infection noise ----
  int<lower=0, upper=1> I_sample; // Stochastic (=1) or deterministic (=0) renewal process?
  int<lower=0, upper=I_sample> I_overdispersion; // whether to model overdispersion via a negative binomial
  array[I_overdispersion ? 2 : 0] real I_xi_prior; // prior on the overdispersion parameter

  // Effective reproduction number ----
  int<lower=0, upper=5> R_model; // 0 for exponential smoothing (ets), 1 for spline smoothing (bs), 2 for soft changepoints (scp)
  int<lower=0, upper=1> R_use_ets;
  int<lower=0, upper=1> R_use_bs;
  int<lower=0, upper=1> R_use_bs2;
  int<lower=0, upper=1> R_use_scp;
  int<lower=0, upper=1> R_use_gp;

  array[2] real R_intercept_prior; // prior for R at start of modeled time series

  // exponential smoothing (ets), random walk is just a special case of ets
  int<lower=0> ets_length; // L + S + D + T - (G+se)
  array[R_use_ets ? 2 : 0] real ets_trend_start_prior;
  array[R_use_ets ? 1 : 0] int<lower=0> ets_diff; // order of differencing
  array[R_use_ets ? 1 : 0] int<lower=0, upper=1> ets_noncentered; // use non-centered parameterization?
  array[R_use_ets ? 2 : 0] real<lower=0> ets_alpha_prior;
  array[R_use_ets ? 2 : 0] real<lower=0> ets_beta_prior;
  array[R_use_ets ? 2 : 0] real<lower=0> ets_phi_prior;

  // basis splines (bs)
  // Sparse bs matrix: columns = bases (bs_ncol), rows = time points (L+S+T-(G+se))
  // Global
  int<lower=0> bs_length; // (L + S + D + T - (G+se) + h)
  int<lower=0> bs_dist; // standard distance between knots
  array[R_use_bs ? 1 : 0] int<lower=1> bs_ncol; // number of B-splines
  array[R_use_bs ? 1 : 0] int<lower=0> bs_n_w; // number of nonzero entries in bs matrix
  vector[R_use_bs ? bs_n_w[1] : 0] bs_w; // nonzero entries in bs matrix
  array[R_use_bs ? bs_n_w[1] : 0] int bs_v; // column indices of bs_w
  array[R_use_bs ? (bs_length + 1) : 0] int bs_u; // row starting indices for bs_w plus padding
  // Local
  int<lower=0> bs2_length; // (L + S + D + T - (G+se) + h)
  int<lower=0> bs2_dist; // standard distance between knots
  array[R_use_bs2 ? 1 : 0] int<lower=1> bs2_ncol; // number of B-splines
  array[R_use_bs2 ? 1 : 0] int<lower=0> bs2_n_w; // number of nonzero entries in bs matrix
  vector[R_use_bs2 ? bs2_n_w[1] : 0] bs2_w; // nonzero entries in bs matrix
  array[R_use_bs2 ? bs2_n_w[1] : 0] int bs2_v; // column indices of bs2_w
  array[R_use_bs2 ? (bs2_length + 1) : 0] int bs2_u; // row starting indices for bs2_w plus padding

  // soft changepoints (scp)
  int<lower=0> scp_length; // L + S + D + T - (G+se)
  array[R_use_scp ? 1 : 0] int<lower=1> scp_n_knots;
  array[R_use_scp ? 1 : 0] int<lower=1> scp_break_dist;
  array[R_use_scp ? 1 : 0] int<lower=1> scp_min_dist;
  array[R_use_scp ? 1 : 0] real<lower=0> scp_skip_tolerance; // tolerance for skipping a changepoint
  array[R_use_scp ? 1 : 0] real<lower=0> scp_skip_tolerance_k; // trength of logistic link for skip tolerance
  array[R_use_scp ? 1 : 0] int<lower=0> scp_length_intercept;
  array[R_use_scp ? 1 : 0] real scp_boltzmann_sharpness; // strength of smooth maximum function for fuzzy OR
  vector[R_use_scp ? scp_break_dist[1] : 0] scp_alpha_base;
  vector[R_use_scp ? scp_break_dist[1] : 0] scp_alpha_adjusted;

  // gaussian process (gp)
  array[R_use_gp ? 1 : 0] int<lower=1> gp_n; // number of time points (L + S + D + T - (G+se) + h - 1)
  array[R_use_gp ? 1 : 0] real<lower=0> gp_matern_nu; // smoothness (Matern kernel)
  array[R_use_gp ? 2 : 0] real<lower=0> gp_sigma_prior; // magnitude (Matern kernel)
  array[R_use_gp ? 2 : 0] real<lower=0> gp_length_prior; // length scale (Matern kernel)
  array[R_use_gp ? 1 : 0] real<lower=0> gp_length_max; // max length scale (Matern kernel)
  array[R_use_gp ? 1 : 0] real<lower=0> gp_c; // boundary factor
  array[R_use_gp ? 1 : 0] int<lower=0> gp_m; // number of basis functions
  array[R_model == 5 ? 1 : 0] real<lower=0> gp_ar_phi; // autoregressive smoothing parameter

  // Change point model for R variability
  array[(R_model == 0 || R_model == 1) ? 2 : 0] real R_sd_baseline_prior; // sd of half-normal prior on baseline R variability
  array[2] real R_sd_change_prior; // shape and rate of lomax prior (exponential with gamma distributed rate) on additive R variability at changepoints
  array[(R_model == 0 || R_model == 1) ? 1 : 0] int<lower=1> R_vari_ncol; // number of B-splines (degree 1) for change points
  array[(R_model == 0 || R_model == 1) ? 1 : 0] int<lower=0> R_vari_n_w; // number of nonzero entries in R_vari matrix
  vector[(R_model == 0 || R_model == 1) ? R_vari_n_w[1] : 0] R_vari_w; // nonzero entries in R_vari matrix
  array[(R_model == 0 || R_model == 1) ? R_vari_n_w[1] : 0] int R_vari_v; // column indices of R_vari_w
  array[(R_model == 0 || R_model == 1) ? L + S + D + T - (G+se) + 1 + h : 0] int R_vari_u; // row starting indices for R_vari_w plus padding

  // Selection matrix for R variability (needed for spline smoothing)
  array[R_model == 1 ? 1 : 0] int<lower=1> R_vari_sel_ncol; // number of columns
  array[R_model == 1 ? 1 : 0] int<lower=0> R_vari_sel_n_w; // number of nonzero entries in R_vari_sel matrix
  vector[R_model == 1 ? R_vari_sel_n_w[1] : 0] R_vari_sel_w; // nonzero entries in R_vari_sel matrix
  array[R_model == 1 ? R_vari_sel_n_w[1]: 0] int R_vari_sel_v; // column indices of R_vari_sel_w
  array[R_model == 1 ? bs_ncol[1] : 0] int R_vari_sel_u; // row starting indices for R_vari_sel_w plus padding
  // for local
  array[R_model == 1 ? 1 : 0] int<lower=1> R_vari_sel_local_ncol; // number of columns
  array[R_model == 1 ? 1 : 0] int<lower=0> R_vari_sel_local_n_w; // number of nonzero entries in R_vari_sel_local matrix
  vector[R_model == 1 ? R_vari_sel_local_n_w[1] : 0] R_vari_sel_local_w; // nonzero entries in R_vari_sel_local matrix
  array[R_model == 1 ? R_vari_sel_local_n_w[1]: 0] int R_vari_sel_local_v; // column indices of R_vari_sel_local_w
  array[R_model == 1 ? bs2_ncol[1] : 0] int R_vari_sel_local_u; // row starting indices for R_vari_sel_local_w plus padding

  // Link function and corresponding hyperparameters ----
  // first element: 0 = inv_softplus, 1 = scaled_logit
  // other elements: hyperparameters for the respective link function
  array[4] real R_link;

  // forecast dampening
  real<lower=0, upper=1> forecast_dampening; // dampening of forecasted R
}
transformed data {
  vector[G] gi_rev = reverse(generation_dist);
  vector[L + 1] inc_rev = reverse(incubation_dist);
  vector[S + 1] shed_rev_log = log(reverse(shedding_dist));
  vector[D + 1] residence_rev_log = log(reverse(residence_dist));
  vector[T+h] flow_log = log(flow);

  real flow_median_log = log(quantile(flow, 0.5));

  // number of averaged technical replicates per date
  real n_averaged_median = quantile(n_averaged, 0.5);
  vector[T+h] n_averaged_all = rep_vector(n_averaged_median, T+h);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the number of replicates of the last one will be used for that date
    n_averaged_all[sample_to_date[measure_to_sample[i]]] = n_averaged[i];
  }

  // number of total partitions per measurement per date
  real total_partitions_median = quantile(dPCR_total_partitions, 0.5);
  vector[T+h] total_partitions_all = rep_vector(total_partitions_median, T+h);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the total partitions of the last one will be used for that date
    total_partitions_all[sample_to_date[measure_to_sample[i]]] = dPCR_total_partitions[i];
  }

  // Upper relevant bound for LOD model
  real conc_drop_prob;
  if (LOD_model == 0) {
    conc_drop_prob = positive_infinity();
  } else {
    real LOD_expected_scale;
    if (LOD_model == 1) {
      LOD_expected_scale = LOD_scale[1];
    } else if (LOD_model == 2) {
      LOD_expected_scale = (
      (total_partitions_observe ? total_partitions_median : (nu_upsilon_b_mu_prior[1] * 1e4)) *
      nu_upsilon_c_prior[1] * 1e-5 *
      n_averaged_median
      );
    }
    conc_drop_prob = -log(LOD_drop_prob)/LOD_expected_scale; // concentrations above this value are irrelevant for LOD model (probability of non-detection is virtually zero)
  }

  int n_zero = num_zeros(measured_concentrations);
  int n_dropLOD = num_zeros(fmax(0, conc_drop_prob - measured_concentrations));
  array[n_measured] int<lower=0> i_all;
  array[n_zero] int<lower=0> i_zero;
  array[n_measured - n_zero] int<lower=0> i_nonzero;
  array[n_measured - n_dropLOD] int<lower=0> i_LOD;
  array[n_measured - n_zero - n_dropLOD] int<lower=0> i_nonzero_LOD;
  int i_z = 0;
  int i_nz = 0;
  int i_lod = 0;
  int i_nzs = 0;
  for (n in 1:n_measured) {
    i_all[n] = n;
    if (measured_concentrations[n] == 0) {
      i_z += 1;
      i_zero[i_z] = n;
    } else {
      i_nz += 1;
      i_nonzero[i_nz] = n;
      if (measured_concentrations[n] < conc_drop_prob) {
        i_nzs += 1;
        i_nonzero_LOD[i_nzs] = n;
      }
    }
    if (measured_concentrations[n] < conc_drop_prob) {
        i_lod += 1;
        i_LOD[i_lod] = n;
    }
  }
  int n_dropzero = (obs_dist == 3 ? 0 : n_zero); // only include zeros for non-truncated normal
  array[n_measured - n_dropzero] int<lower=0> i_include;
  if (n_dropzero > 0) {
    i_include = i_nonzero;
  } else {
    i_include = i_all;
  }

  matrix[R_use_gp ? gp_n[1] : 0, R_use_gp ? gp_m[1] : 0] gp_PHI;
  array[R_use_gp ? 1 : 0] real<lower=0> gp_L; // boundary condition
  int gp_params_fixed = R_use_gp && gp_sigma_prior[2] == 0 && gp_length_prior[2] == 0;
  vector[gp_params_fixed ? gp_m[1] : 0] diagSPD_fixed;
  if (R_use_gp) {
    // maximum absolute value of input space (mean-centered time points)
    real gp_S = (gp_n[1] - 1) / 2.0;
    // mean-centered time points in interval [-S, S]
    vector[gp_n[1]] x = linspaced_vector(gp_n[1], -gp_S, gp_S);
    // boundary condition = boundary factor * S
    gp_L[1] = gp_c[1] * gp_S;
    // eigenvalue matrix of Laplacian operator in the interval [-L, L]
    gp_PHI = gp_eigenvalue_matrix(gp_n[1], gp_m[1], gp_L[1], x);
    if (gp_params_fixed) {
      // if parameters are fixed, precompute diagonal of SPD matrix
      diagSPD_fixed = diagSPD_Matern(
        gp_matern_nu[1], // smoothness
        gp_sigma_prior[1], // magnitude
        gp_length_prior[1], // length scale
        gp_L[1], // boundary condition
        gp_m[1] // number of basis functions
      );
    }
  }

  if (R_link[1] < 0 || R_link[1] > 1) {
    reject("Link function must be one of inv_softplus (0) or scaled_logit (1)");
  }
}
parameters {
  // Effective reproduction number parameters ----
  real R_intercept; // starting value of R

  // random walk / exponential smoothing (ets)
  array[R_use_ets ? 1 : 0] real ets_trend_start; // starting value of the trend
  vector[R_use_ets ? ets_length - 1 : 0] ets_noise; // additive errors
  array[R_use_ets && ets_alpha_prior[2] > 0 ? 1 : 0] real<lower=0, upper=1> ets_alpha; // smoothing parameter for the level
  array[R_use_ets && ets_beta_prior[2] > 0 ? 1 : 0] real<lower=0, upper=1> ets_beta; // smoothing parameter for the trend
  array[R_use_ets && ets_phi_prior[2] > 0 ? 1 : 0] real<lower=0, upper=1> ets_phi; // dampening parameter of the trend

  // basis splines (bs)
  vector[R_use_bs ? (R_model == 4 ? bs_ncol[1] - 2 : bs_ncol[1] - 1) : 0] bs_coeff_noise_raw; // additive errors (non-centered)
  vector[R_use_bs2 ? (bs2_ncol[1] - 1) : 0] bs2_coeff_noise_raw; // additive errors (non-centered)

  // smooth derivative
  vector<lower=0>[R_model == 4 ? bs_ncol[1]-2 : 0] bs_coeff_noise_lomax;

  // soft changepoints (scp)
  array[R_use_scp ? scp_n_knots[1] : 1] vector[R_use_scp ? (scp_break_dist[1]-1) : 1] scp_break_delays_raw;
  vector[R_use_scp ? scp_n_knots[1] : 0] scp_noise; // additive errors
  vector<lower=0>[R_use_scp ? scp_n_knots[1] : 0] scp_sd; // R variability for soft changepoint model

  // Gaussian process (gp)
  vector<lower=0>[R_model == 5 ? gp_m[1] : 0] gp_noise_raw; // non-centered noise for basis functions
  array[R_use_gp && gp_sigma_prior[2] > 0 ? 1 : 0] real<lower=0> gp_sigma; // magnitude
  array[R_use_gp && gp_length_prior[2] > 0 ? 1 : 0] real<lower=0, upper=(R_use_gp ? gp_length_max[1] : 0)> gp_length; // length scale

  // Change point model for Rt variability
  array[(R_model == 0 || R_model == 1) ? (R_sd_baseline_prior[2] > 0 ? 1 : 0) : 0] real<lower=0> R_sd_baseline; // baseline R variability
  vector<lower=0>[(R_model == 0 || R_model == 1) ? (R_vari_ncol[1] > 2 ? R_vari_ncol[1] - 2 : R_vari_ncol[1]) : 0] R_sd_changepoints; // additive R variability at changepoints

  // seeding
  real iota_log_seed_intercept;
  array[seeding_model > 0 ? 1 : 0] real<lower=0> iota_log_seed_sd;
  vector<multiplier=(seeding_model > 0 ? iota_log_seed_sd[1] : 1)>[seeding_model > 0 ? (G+se) - 1 : 0] iota_log_ar_noise;

  // realized infections
  array[I_overdispersion && (I_xi_prior[2] > 0) ? 1 : 0] real<lower=0> I_xi; // positive to ensure identifiability
  vector<lower=0>[I_sample ? L + S + D + T : 0] I; // realized number of infections

  // paramaters of shedding load distribution
  simplex[shedding_dist_n > 1 ? shedding_dist_n : 1] shedding_dist_weights;
  vector<lower=0>[shedding_dist_type > 0 ? shedding_dist_n : 0] shedding_dist_mean;
  vector<lower=0>[shedding_dist_type > 1 ? shedding_dist_n : 0] shedding_dist_cv;

  // individual-level shedding load variation
  array[load_vari && nu_zeta_prior[2] > 0 ? 1 : 0] real<lower=0> nu_zeta; // coefficient of variation of individual-level load
  vector[load_vari ? n_zeta_exact : 0] zeta_log_exact; // realized shedding load
  vector[load_vari ? n_zeta_normal_approx : 0] zeta_raw_approx; // realized shedding load (non-centered)

  // sample date effects
  vector[K] eta;

  // Coefficient of variation of likelihood for measurements
  array[pr_noise ? 1 : 0] real<lower=0> nu_psi; // pre-replicaton coefficient of variation
  vector<multiplier = (pr_noise ? nu_psi[1] : 1)>[pr_noise ? n_samples : 0] psi; // realized concentration before replication stage
  real<lower=0> nu_upsilon_a;
  array[(cv_type == 1) && (nu_upsilon_b_mu_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_mu;
  array[(cv_type == 1) && total_partitions_observe!=1 && (nu_upsilon_b_cv_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_cv;
  vector[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b_noise_raw;
  array[(cv_type == 1) && nu_upsilon_c_prior[2] > 0 ? 1 : 0] real<lower=0> nu_upsilon_c;

  vector<lower=0>[outliers ? D + T : 0] epsilon; // external noise at low concentrations;
}
transformed parameters {
  // Effective reproduction number parameters ----
  vector[L + S + D + T - G] R;
  // random walk / exponential smoothing (ets)
  vector<lower=0>[R_use_ets ? ets_length - 1 : 0] ets_sd; // standard deviation of additive errors in R ets model
  // soft changepoints (scp)
  vector[R_use_scp ? scp_n_knots[1] : 0] scp_knot_values;
  vector[R_model == 3 ? scp_length : 0] scp_values;
  array[R_use_scp ? scp_n_knots[1] : 1] simplex[R_use_scp ? scp_break_dist[1] : 1] scp_break_delays;
  if (1-R_use_scp) {
    scp_break_delays[1] = rep_vector(1, 1); // if no soft changepoints, set to 1
  }

  // parameters for specific R models
  vector[R_model == 1 || R_model == 5 ? h : 0] R_forecast_model; // spline-based forecast of R
  vector<lower=0>[R_model == 1 ? bs_length : 0] bs_coeff_ar_sd; // sd for random walk on log bs coeffs

  // other time series parameters ----
  vector[L + S + D + T] iota; // expected number of infections
  vector[S + D + T] lambda; // expected number of shedding onsets
  vector[D + T] omega_log; // log expected daily loads in catchment
  vector[T] pi_log; // log expected daily loads at sampling site
  vector[T] kappa_log; // log expected daily concentrations
  vector[n_samples] rho_log; // log expected concentrations in (composite) samples

  // shedding-related parameters ----
  vector[shedding_dist_type > 0 ? S + 1 : 0] shed_rev_log_sample; // shedding load distribution
  vector[load_vari ? S + D + T : 0] zeta_log; // realized shedding load

  // dPCR model ----
  vector<lower=0>[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b; // total partitions per measurement
  array[LOD_model > 0 ? 1 : 0] vector<lower=0>[n_measured] LOD_hurdle_scale;

  // Reproduction number ----
  if (R_model == 0) {
    ets_sd = csr_matrix_times_vector(
      L + S + D + T - (G+se) + h, R_vari_ncol[1], R_vari_w,
      R_vari_v, R_vari_u,
      param_or_fixed(R_sd_baseline, R_sd_baseline_prior) +
      (R_vari_ncol[1] > 2 ? append_row2([0]', R_sd_changepoints, [0]') : R_sd_changepoints)
      )[2:(L + S + D + T - (G+se))];
    // Innovations state space process implementing exponential smoothing
    R[(se+1):(L + S + D + T - G)] = apply_link(holt_damped_process(
      [R_intercept, ets_trend_start[1]]',
      param_or_fixed(ets_alpha, ets_alpha_prior),
      param_or_fixed(ets_beta, ets_beta_prior),
      param_or_fixed(ets_phi, ets_phi_prior),
      ets_noise .* (ets_noncentered[1] ? ets_sd : rep_vector(1, ets_length - 1)),
      ets_diff[1]
    ), R_link);
  } else if (R_model == 1) {
    // Basis spline smoothing
    // global
    vector[bs_length] R_global;
    bs_coeff_ar_sd = csr_matrix_times_vector(
      bs_length, R_vari_ncol[1], R_vari_w,
      R_vari_v, R_vari_u,
      (R_vari_ncol[1] > 2 ? append_row2(
        [0]', R_sd_changepoints, [0]'
        ) : R_sd_changepoints)
      );
    vector[bs_ncol[1] - 1] bs_coeff_noise = bs_coeff_noise_raw .*
    sqrt(csr_matrix_times_vector(
      bs_ncol[1] - 1, R_vari_sel_ncol[1], R_vari_sel_w,
      R_vari_sel_v, R_vari_sel_u,
      bs_coeff_ar_sd^2 // we add together variances --> square sd
      ));
    vector[bs_ncol[1]] bs_coeff = random_walk([0]', bs_coeff_noise, 0); // Basis spline coefficients
    bs_coeff[(bs_ncol[1]-3):bs_ncol[1]] = rep_vector(bs_coeff[bs_ncol[1]-4],4); // fix forecast at last value
    R_global = csr_matrix_times_vector(
      bs_length, bs_ncol[1], bs_w, bs_v, bs_u, bs_coeff
      );
    // local
    vector[bs2_length] R_local;
    vector[bs2_ncol[1] - 1] bs2_coeff_noise = bs2_coeff_noise_raw .*
    sqrt(csr_matrix_times_vector(
      bs2_ncol[1] - 1, R_vari_sel_local_ncol[1], R_vari_sel_local_w,
      R_vari_sel_local_v, R_vari_sel_local_u,
      rep_vector(param_or_fixed(R_sd_baseline, R_sd_baseline_prior)^2, bs2_length) // we add together variances --> square sd
      ));
    vector[bs2_ncol[1]] bs2_coeff = random_walk([0]', bs2_coeff_noise, 0); // Basis spline coefficients
    R_local = csr_matrix_times_vector(
      bs2_length, bs2_ncol[1], bs2_w, bs2_v, bs2_u, bs2_coeff
      );
    // intercept + global + local
    vector[L + S + D + T - (G+se) + h] R_all = apply_link(
      R_intercept + R_global + R_local, R_link
      );
    R[(se+1):(L + S + D + T - G)] = R_all[1:(L + S + D + T - (G+se))];
    if (h>0) {
      R_forecast_model = R_all[(L + S + D + T - (G+se) + 1):(L + S + D + T - (G+se) + h)];
    }
  } else if (R_model == 2) {
    for (i in 1:scp_n_knots[1]) {
      scp_break_delays[i] = softmax(append_row(scp_break_delays_raw[i], 0));
    }
    scp_knot_values = random_walk(
      [R_intercept]', scp_noise .* scp_sd, 0
      )[2:(scp_n_knots[1]+1)]; // do not keep intercept here
    R[(se+1):(L + S + D + T - G)] = apply_link(soft_changepoint(
      R_intercept, scp_knot_values, scp_length_intercept[1], scp_length,
      scp_break_dist[1], scp_min_dist[1], scp_break_delays,
      scp_boltzmann_sharpness[1], scp_skip_tolerance[1], scp_skip_tolerance_k[1]
      ), R_link);
  } else if (R_model == 3) {
    for (i in 1:scp_n_knots[1]) {
      scp_break_delays[i] = softmax(append_row(scp_break_delays_raw[i], 0));
    }
    scp_knot_values = scp_noise .* scp_sd;
    scp_values = cumulative_sum(soft_changepoint(
      0, scp_knot_values, scp_length_intercept[1], scp_length,
      scp_break_dist[1], scp_min_dist[1], scp_break_delays,
      scp_boltzmann_sharpness[1], scp_skip_tolerance[1], scp_skip_tolerance_k[1]
      ));
    real last_diff = scp_knot_values[scp_n_knots[1]];
    vector[bs_ncol[1]] bs_coeff = append_row2(
      [scp_values[1]]',
      scp_values,
      [
      scp_values[scp_length] + last_diff * (1 - 1/3.0),
      scp_values[scp_length] + last_diff
      ]'
      );
    R[(se+1):(L + S + D + T - G)] = apply_link(
      R_intercept + csr_matrix_times_vector(
        bs_length, bs_ncol[1], bs_w, bs_v, bs_u, bs_coeff
      ), R_link);
  } else if (R_model == 4) {
    vector[bs_ncol[1]] bs_coeff = append_row2(
      0, (sqrt(bs_coeff_noise_lomax) .* bs_coeff_noise_raw), 0
      );
    vector[(L + S + D + T - G) - se] R_spline = csr_matrix_times_vector(
      bs_length, bs_ncol[1], bs_w, bs_v, bs_u, bs_coeff
      );
    R[(se+1):(L + S + D + T - G)] = apply_link(
      R_intercept + cumulative_sum(R_spline), R_link
      );
  } else if (R_model == 5) {
    vector[gp_m[1]] diagSPD;
    if (gp_params_fixed) {
      diagSPD = diagSPD_fixed;
    } else {
      diagSPD = diagSPD_Matern(
        gp_matern_nu[1], // smoothness
        param_or_fixed(gp_sigma, gp_sigma_prior), // magnitude
        param_or_fixed(gp_length, gp_length_prior), // length scale
        gp_L[1], // boundary condition
        gp_m[1] // number of basis functions
        );
    }
    vector[L + S + D + T - (G+se) + h] R_gp = ar1_process(
      R_intercept - 1, gp_ar_phi[1], gp_PHI * (diagSPD .* gp_noise_raw)
      );
    vector[L + S + D + T - (G+se) + h] R_all = apply_link(1 + R_gp, R_link);
    R[(se+1):(L + S + D + T - G)] = R_all[1:(L + S + D + T - (G+se))];
    if (h>0) {
      R_forecast_model = R_all[(L + S + D + T - (G+se) + 1):(L + S + D + T - (G+se) + h)];
    }
  }

  // seeding
  if (seeding_model == 0) {
    iota[1:(G+se)] = exp(
      rep_vector(iota_log_seed_intercept, (G+se))
      );
  } else if (seeding_model == 1) {
    iota[1:(G+se)] = exp(
      random_walk([iota_log_seed_intercept]', iota_log_ar_noise, 0)
      );
  } else if (seeding_model == 2) {
    real seeding_last_r = get_growth_rate(R[se+1], generation_dist);
    vector[G + se] seeding_r = reverse(
      random_walk([seeding_last_r]', iota_log_ar_noise, 0)
      ); // backward-in-time random walk on growth rate
    iota[1:(G+se)] = exp(
      random_walk([iota_log_seed_intercept]', seeding_r[1:(G+se-1)], 0)
      );
  }

  // compute Rt for extended seeding phase
  if (se > 0) {
    vector[L + S + D + T - G] infness;
    infness = infectiousness(L + S + D + T - G, G, 0, gi_rev, iota);
    R[1:se] = iota[(G+1):(G+se)] ./ infness[1:se];
  }

  // renewal process
  if (I_sample) {
    iota[(G + se + 1) : (L + S + D + T)] = renewal_process_stochastic(
      (L + S + D + T - (G+se)), R[(se+1):(L + S + D + T - G)], G, se, gi_rev, I);
  } else {
    iota[(G + se + 1) : (L + S + D + T)] = renewal_process_deterministic(
      (L + S + D + T - (G+se)), R[(se+1):(L + S + D + T - G)], G, se, gi_rev, iota);
  }

  // convolution from infections to shedding onsets (expected)
  lambda = convolve(inc_rev, I_sample ? I : iota)[(L + 1) : (L + S + D + T)];

  // sampled shedding distribution
  if (shedding_dist_type > 0) {
    shed_rev_log_sample = reverse(
      discretise_dist_log(
        shedding_dist_type,
        shedding_dist_n > 1 ? dot_product(shedding_dist_weights, shedding_dist_mean) : shedding_dist_mean[1],
        shedding_dist_type > 1 ?  (shedding_dist_n > 1 ? dot_product(shedding_dist_weights, shedding_dist_cv) : shedding_dist_cv[1]) : 0,
        S
      )
    );
  }

  // calculation of total loads shed each day (expected)
  if (load_vari) {
    // Shedding load variation
    zeta_log[zeta_exact] = zeta_log_exact;
    zeta_log[zeta_normal_approx] = gamma_sum_log_approx(
      param_or_fixed(nu_zeta, nu_zeta_prior),
      lambda[zeta_normal_approx],
      zeta_raw_approx
    );
    omega_log = log_convolve(
        shedding_dist_type == 0 ? shed_rev_log : shed_rev_log_sample, // shedding load distribution
        log(load_mean) + zeta_log // total load shed
        )[(S + 1) : (S + D + T)];
  } else {
    omega_log = log_convolve(
        shedding_dist_type == 0 ? shed_rev_log : shed_rev_log_sample, // shedding load distribution
        log(load_mean) + log(lambda) // total load shed
        )[(S + 1) : (S + D + T)];
  }

  // calculation of total loads at sampling site (expected)
  if (D>0) {
    pi_log = log_convolve(
      residence_rev_log, // residence time distribution
      omega_log
      )[(D + 1) : (D + T)];
  } else {
    pi_log = omega_log;
  }

  // calculation of concentrations at measurement site by day (expected)
  // --> adjusted for flow and for date of sample effects
  if (K > 0) {
    kappa_log = pi_log - flow_log[1:T] + X[1:T] * eta;
  } else {
    kappa_log = pi_log - flow_log[1:T];
  }

  if (outliers) {
    kappa_log = log_sum_exp(kappa_log, log(load_mean) + log(epsilon) - flow_median_log); // additive outlier component
  }

  // concentrations in (composite) samples
  if (w > 1) { // multi-day composite samples
    for (i in 1 : n_samples) {
        rho_log[i] = log_sum_exp(
          kappa_log[(sample_to_date[i] - w + 1) : sample_to_date[i]]
          ) - log(w);
    }
  } else { // individual day samples
     for (i in 1 : n_samples) {
        rho_log[i] = kappa_log[sample_to_date[i]];
     }
  }

  if ((cv_type == 1) && total_partitions_observe!=1) {
    nu_upsilon_b = total_partitions_noncentered(
      param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
      param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
      nu_upsilon_b_noise_raw
    );
  }

  // scale for LOD hurdle model
  if (LOD_model == 1) {
    LOD_hurdle_scale[1] = rep_vector(LOD_scale[1], n_measured);
  } else if (LOD_model == 2) {
    LOD_hurdle_scale[1] = (
    n_averaged .*
    (total_partitions_observe ? dPCR_total_partitions : nu_upsilon_b * 1e4) *
     param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
    );
  }
}
model {
  // Priors
  R_intercept ~ normal(R_intercept_prior[1], R_intercept_prior[2]); // prior for Rt at start of time series
  if (R_model == 0 || R_model == 1) {
    target += normal_prior_lpdf(R_sd_baseline | R_sd_baseline_prior, 0); // truncated normal
    target += pareto_type_2_lpdf(
      R_sd_changepoints | 0, R_sd_change_prior[2], R_sd_change_prior[1]
      ); // Lomax distribution / exponential-gamma distribution
  } else if (R_model == 2 || R_model == 3) {
    target += pareto_type_2_lpdf(
      scp_sd | 0, R_sd_change_prior[2], R_sd_change_prior[1]
    );
    scp_noise ~ std_normal();
  } else if (R_model == 4) {
    target += pareto_type_2_lpdf(
      bs_coeff_noise_lomax | 0, R_sd_change_prior[2], R_sd_change_prior[1]
    );
  } else if (R_model == 5) {
    target += normal_prior_lpdf(gp_sigma | gp_sigma_prior, 0); // truncated normal
    target += normal_prior_lb_ub_lpdf(gp_length | gp_length_prior, 0, gp_length_max[1]); // truncated normal
    gp_noise_raw ~ std_normal(); // Gaussian noise for GP model
  }

  if (R_use_ets) {
    // Innovations state space model
    target += ets_coefficient_priors_lp(
      ets_alpha, ets_alpha_prior,
      ets_beta, ets_beta_prior,
      ets_phi, ets_phi_prior
      );
    ets_trend_start ~ normal(ets_trend_start_prior[1], ets_trend_start_prior[2]); // starting prior for trend
    if (ets_noncentered[1]) {
      ets_noise ~ std_normal(); // Gaussian noise
    } else {
      ets_noise ~ normal(0, ets_sd); // Gaussian noise
    }
  }
  if (R_use_bs) {
    bs_coeff_noise_raw ~ std_normal(); // Gaussian noise
  }
  if (R_use_bs2) {
    bs2_coeff_noise_raw ~ std_normal(); // Gaussian noise
  }
  if (R_use_scp) {
    target += dirichlet_reduced_lpdf(
      scp_break_delays_raw[1] | scp_alpha_base
    );
    for (i in 2:scp_n_knots[1]) {
      target += dirichlet_reduced_lpdf(
        scp_break_delays_raw[i] | scp_alpha_adjusted
      );
    }
  }

  // Seeding
  iota_log_seed_intercept ~ normal(iota_log_seed_intercept_prior[1], iota_log_seed_intercept_prior[2]);
  if (seeding_model > 0) {
    iota_log_seed_sd[1] ~ normal(iota_log_seed_sd_prior[1], iota_log_seed_sd_prior[2]) T[0, ]; // truncated normal
    iota_log_ar_noise ~ normal(0, iota_log_seed_sd[1]); // Gaussian noise
  }

  // Sampling of infections
  if (I_sample) {
    if (I_overdispersion) {
      target += normal_prior_lpdf(I_xi | I_xi_prior, 0); // truncated normal
      I[1 : (L + S + D + T)] ~ normal(iota, iota .* soft_upper(sqrt(iota .* (1 + iota * (param_or_fixed(I_xi, I_xi_prior) ^ 2)))./iota, 0.5, 10)) T[0, ]; // approximates negative binomial
    } else {
      I[1 : (L + S + D + T)] ~ normal(iota, iota .* soft_upper(sqrt(iota)./iota, 0.5, 10)) T[0, ]; // approximates Poisson
    }
  }

  // Parameters of shedding load distribution
  if (shedding_dist_type > 0) {
    if (shedding_dist_n > 1) {
    shedding_dist_weights ~ dirichlet(shedding_dist_weights_prior);
    }
    for (i in 1:shedding_dist_n) {
      shedding_dist_mean[i] ~ normal(shedding_dist_mean_prior[i,1], shedding_dist_mean_prior[i,2]) T[0, ]; // truncated normal
      if (shedding_dist_type > 1) {
        shedding_dist_cv[i] ~ normal(shedding_dist_cv_prior[i,1], shedding_dist_cv_prior[i,2]) T[0, ]; // truncated normal
      }
    }
  }

  // Prior on individual-level shedding load variation
  if (load_vari) {
    target += normal_prior_lpdf(nu_zeta | nu_zeta_prior, 0); // truncated normal
    target += gamma3_sum_log_lpdf(zeta_log_exact | 1, param_or_fixed(nu_zeta, nu_zeta_prior), lambda[zeta_exact]); // gamma sum (exact)
    zeta_raw_approx ~ std_normal(); // non-centered noise for gamma sum (approx)
  }

  // Prior for additive outlier component
  if (outliers) {
    target += gev_lpdf(epsilon | epsilon_prior[1], epsilon_prior[2], epsilon_prior[3]); // generalized extreme value distribution
  }

  // Prior on sample date effects
  if (K > 0) {
    eta ~ normal(eta_prior[1], eta_prior[2]);
  }

  // Prior on cv of likelihood for measurements
  nu_upsilon_a ~ normal(nu_upsilon_a_prior[1], nu_upsilon_a_prior[2]) T[0, ]; // truncated normal
  if (cv_type == 1) {
    target += normal_prior_lpdf(nu_upsilon_b_mu | nu_upsilon_b_mu_prior, 0); // truncated normal
    target += normal_prior_lpdf(nu_upsilon_c | nu_upsilon_c_prior, 0); // truncated normal
    if (total_partitions_observe != 1) {
      target += normal_prior_lpdf(nu_upsilon_b_cv | nu_upsilon_b_cv_prior, 0); // truncated normal
      nu_upsilon_b_noise_raw ~ std_normal();
    }
  }

  // Likelihood
  {
    // accounting for pre-replicate noise
    vector[n_measured] concentration_log;
    if (pr_noise) {
      // Prior on cv of pre-replication concentrations
      nu_psi[1] ~ normal(nu_psi_prior[1], nu_psi_prior[2]) T[0, ]; // truncated normal
      // lognormal distribution modeled as normal on the log scale
      // with mean = rho (mean_log = rho_log) and cv = nu_psi
      target += lognormal_log_lpdf(psi | rho_log, nu_psi[1]);
      concentration_log = psi[measure_to_sample];
    } else {
      concentration_log = rho_log[measure_to_sample];
    }

    // concentration on unit scale
    vector[n_measured] concentration = exp(concentration_log);
    vector[n_measured] p_zero_log = rep_vector(negative_infinity(), n_measured);
    vector[n_measured] p_zero = rep_vector(0, n_measured);

    // limit of detection
    if (LOD_model > 0) {
      p_zero_log[i_LOD] = log_hurdle_exponential(
        concentration[i_LOD],
        LOD_hurdle_scale[1][i_LOD], // LOD scale (c * m * n)
        cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
        cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
        );
      target += sum(p_zero_log[i_zero]); // below-LOD probabilities
      target += sum(log1m_exp(p_zero_log[i_nonzero_LOD])); // above-LOD probabilities
      p_zero = exp(p_zero_log);
    }

    // CV of each observation as a function of concentration
    vector[n_measured] cv;
    if (cv_type == 0) { // constant cv
      cv[i_include] = rep_vector(nu_upsilon_a, n_measured - n_zero);
    } else if (cv_type == 1) { // dPCR
      cv[i_include] = cv_dPCR_pre(
        concentration[i_include], // lambda (concentration)
        nu_upsilon_a, // nu_pre (pre-PCR CV)
        (total_partitions_observe ? dPCR_total_partitions[i_include] : nu_upsilon_b[i_include] * 1e4), // m (number of partitions)
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
        n_averaged[i_include], // n (number of averaged replicates)
        cv_pre_type[1], // Type of pre-PCR CV
        cv_pre_approx_taylor[1] // Should taylor approximation be used?
        );
    } else if (cv_type == 2) { // constant variance
      cv[i_include] = (
        nu_upsilon_a * mean(measured_concentrations[i_include]) /
        concentration[i_include]
        );
    }

    vector[n_measured] mean_conditional;
    vector[n_measured] cv_conditional;
    mean_conditional[i_include] = concentration[i_include] ./ (1-p_zero[i_include]);
    cv_conditional[i_include] = sqrt(cv[i_include]^2 .* (1-p_zero[i_include]) - p_zero[i_include]);

    // measurements
    if (obs_dist == 0) {
      target += gamma3_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include] // coefficient of variation
      );
    } else if (obs_dist == 1) {
      target += lognormal5_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include] // coefficient of variation
      );
    } else if (obs_dist == 2) {
      target += normal2_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include], // coefficient of variation,
        0 // truncate at zero
      );
    } else if (obs_dist == 3) {
      target += normal2_lpdf(
        measured_concentrations[i_include] |
        mean_conditional[i_include], // expectation
        cv_conditional[i_include] // coefficient of variation
      );
    } else {
      reject("Distribution not supported.");
    }
  }
}
generated quantities {
  // predicted measurements and forecasts
  // note that we here assume the same measurement variance as from composite samples,
  // which may be smaller than that of hypothetical daily measurements
  vector[T] predicted_concentration;
  vector[T] predicted_concentration_norm; // flow-normalized
  vector[h] predicted_concentration_forecast;
  vector[h] predicted_concentration_forecast_norm; // flow-normalized
  vector[h] R_forecast;
  vector[h] iota_forecast;
  vector[h] I_forecast;
  vector[h] lambda_forecast;
  vector[h] omega_log_forecast;
  vector[h] pi_log_forecast;
  vector[h] kappa_log_forecast;
  {
    vector[T+h] concentration_log;
    vector[T+h] concentration;
    vector[T+h] cv_all;
    vector[T+h] isnonzero; // will be a vector of 0s and 1s

    // Prediction for days until present
    if (pr_noise) {
      concentration_log[1:T] = to_vector(lognormal_log_rng(kappa_log, nu_psi[1]));
    } else {
      concentration_log[1:T] = kappa_log;
    }

    // Forecasting for days beyond present
    if (h>0) {
      // Forecasting of R
      if (R_model == 0) {
        // Innovations state space process implementing exponential smoothing
        vector[h] ets_sd_forecast = csr_matrix_times_vector(
          L + S + D + T - (G+se) + h, R_vari_ncol[1], R_vari_w,
          R_vari_v, R_vari_u,
          param_or_fixed(R_sd_baseline, R_sd_baseline_prior) +
          (R_vari_ncol[1] > 2 ? append_row2([0]', R_sd_changepoints, [0]') : R_sd_changepoints)
        )[((L + S + D + T - (G+se)) + 1):((L + S + D + T - (G+se)) + h)];
        R_forecast = apply_link(holt_damped_process(
          [R_intercept, ets_trend_start[1]]',
          param_or_fixed(ets_alpha, ets_alpha_prior),
          param_or_fixed(ets_beta, ets_beta_prior),
          param_or_fixed(ets_phi, ets_phi_prior),
          append_row(
            ets_noise .* (ets_noncentered[1] ? ets_sd : rep_vector(1, L + S + D + T - (G+se) - 1)),
            to_vector(normal_rng(0, ets_sd_forecast))
            ),
          ets_diff[1]
        ), R_link)[((L + S + D + T - (G+se)) + 1):((L + S + D + T - (G+se)) + h)];
      } else if (R_model == 1) {
        R_forecast = R_forecast_model;
      } else if (R_model == 2) {
        R_forecast = rep_vector(R[L + S + D + T - G], h);
      } else if (R_model == 3) {
        real last_R = R[L + S + D + T - G];
        real slope = scp_knot_values[scp_n_knots[1]] / (1.0*bs_dist);
        R_forecast = last_R + cumulative_sum(rep_vector(slope, h));
      } else if (R_model == 4) {
        R_forecast = rep_vector(R[L + S + D + T - G], h);
      } else if (R_model == 5) {
        R_forecast = R_forecast_model;
      }

      R_forecast = dampen_trend(
        R[L + S + D + T - G], R_forecast, forecast_dampening
        );

      // Forecasting of infections
      if (I_sample) {
        array[2] vector[G+h] forecast_tmp;
        forecast_tmp = renewal_process_stochastic_sim_rng(
          h, R_forecast, G, 0, gi_rev,
          append_row(iota[((L + S + D + T)+1-G):(L + S + D + T)], rep_vector(0, h)),
          I_overdispersion ? param_or_fixed(I_xi, I_xi_prior) : 0
          );
        iota_forecast = forecast_tmp[1][(G+1):(G+h)];
        I_forecast = forecast_tmp[2][(G+1):(G+h)];
      } else {
        iota_forecast = renewal_process_deterministic(
          h, R_forecast, G, 0, gi_rev, append_row(iota[((L + S + D + T)+1-G):(L + S + D + T)], rep_vector(0, h))
          );
        I_forecast = iota_forecast;
      }

      // Forecasting of symptom onsets
      lambda_forecast = convolve(
        inc_rev, append_row((I_sample ? I : iota)[((L + S + D + T)+1-L):(L + S + D + T)], I_forecast)
        )[(L+1):(L+h)];

      // Forecasting of total loads
      if (load_vari) {
        vector[h] zeta_log_forecast = log(to_vector(gamma_rng(
          lambda_forecast/param_or_fixed(nu_zeta, nu_zeta_prior)^2,
          1/param_or_fixed(nu_zeta, nu_zeta_prior)^2
          )));
        omega_log_forecast = log_convolve(
            shedding_dist_type == 0 ? shed_rev_log : shed_rev_log_sample, // shedding load distribution
            log(load_mean) + append_row(zeta_log[((S + D + T)+1-S):(S + D + T)], zeta_log_forecast) // total load shed
            )[(S+1):(S+h)];
      } else {
        omega_log_forecast = log_convolve(
            shedding_dist_type == 0 ? shed_rev_log : shed_rev_log_sample, // shedding load distribution
            log(load_mean) + log(append_row(lambda[((S + D + T)+1-S):(S + D + T)], lambda_forecast)) // total load shed
            )[(S+1):(S+h)];
      }

      // Forecasting of total loads at sampling site (expected)
      if (D>0) {
        pi_log_forecast = log_convolve(
          residence_rev_log, // residence time distribution
          append_row(omega_log[((D + T)+1-D):(D + T)], omega_log_forecast)
          )[(D+1):(D+h)];
      } else {
        pi_log_forecast = omega_log_forecast;
      }

      // Forecasting of concentrations at measurement site by day (expected)
      // --> adjusted for flow and for date of sample effects
      if (K > 0) {
        kappa_log_forecast = pi_log_forecast - flow_log[(T+1):(T+h)] + X[(T+1):(T+h)] * eta;
      } else {
        kappa_log_forecast = pi_log_forecast - flow_log[(T+1):(T+h)];
      }

      // Forecasting of pre-replication concentrations
      if (pr_noise) {
        concentration_log[(T+1):(T+h)] = to_vector(lognormal_log_rng(kappa_log_forecast, nu_psi[1]));
      } else {
        concentration_log[(T+1):(T+h)] = kappa_log_forecast;
      }
    }

    concentration = exp(concentration_log);

    vector[T+h] nu_upsilon_b_all;
    if (cv_type == 1) {
      if (total_partitions_observe == 1) {
        nu_upsilon_b_all = total_partitions_all / 1e4;
      } else {
        nu_upsilon_b_all = total_partitions_noncentered(
          param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
          param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
          std_normal_n_rng(T+h)
          );
        for (i in 1:n_measured) {
          nu_upsilon_b_all[sample_to_date[measure_to_sample[i]]] = nu_upsilon_b[i];
        }
      }
    }

    vector[T+h] p_zero_all;
    if (LOD_model > 0) {
      vector[T+h] LOD_hurdle_scale_all;
      // scale for LOD hurdle model
      if (LOD_model == 1) {
        LOD_hurdle_scale_all = rep_vector(LOD_scale[1], T+h);
      } else if (LOD_model == 2) {
        LOD_hurdle_scale_all = (
        n_averaged_all .*
        nu_upsilon_b_all * 1e4 *
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
        );
      }
      p_zero_all = exp(log_hurdle_exponential(
          concentration, // lambda (concentration)
          LOD_hurdle_scale_all,
          cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
          cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
          ));
      p_zero_all = trim_or_reject_ub(
        p_zero_all,
        1-1e-5, // trim to almost 1
        1.01 // throw error when significantly above 1
      );
      isnonzero = to_vector(bernoulli_rng(1-p_zero_all));
    } else {
      p_zero_all = rep_vector(0, T+h);
      isnonzero = rep_vector(1, T+h);
    }

    if (cv_type == 0) {
      cv_all = rep_vector(nu_upsilon_a, T+h);
    } else if (cv_type == 1) {
      cv_all = cv_dPCR_pre(
        concentration, // lambda (concentration)
        nu_upsilon_a, // nu_pre (pre-PCR CV)
        nu_upsilon_b_all * 1e4, // m (number of partitions)
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
        n_averaged_all, // n (number of averaged replicates) for all dates
        cv_pre_type[1], // Type of pre-PCR CV
        cv_pre_approx_taylor[1] // Should taylor approximation be used?
        );
    } else if (cv_type == 2) {
      cv_all = (
        nu_upsilon_a * mean(measured_concentrations[i_nonzero]) /
        concentration
        );
    }

    vector[T+h] mean_conditional_all = concentration ./ (1-p_zero_all);
    vector[T+h] cv_conditional_all = sqrt(cv_all^2 .* (1-p_zero_all) - p_zero_all);

    // correct potentially slightly negative approximations
    cv_conditional_all = trim_or_reject_lb(
      cv_conditional_all,
      1e-5, // trim to almost zero
      -0.01 // throw error when significantly below zero
    );

    vector[T+h] meas_conc;
    if (obs_dist == 0) {
      meas_conc = gamma3_rng(mean_conditional_all, cv_conditional_all);
    } else if (obs_dist == 1) {
      meas_conc = lognormal5_rng(mean_conditional_all, cv_conditional_all);
    } else if (obs_dist == 2) {
      meas_conc = normal2_rng(mean_conditional_all, cv_conditional_all, 0); // truncated at zero
    } else if (obs_dist == 3) {
      meas_conc = normal2_rng(mean_conditional_all, cv_conditional_all);
    } else {
      reject("Distribution not supported.");
    }
    predicted_concentration = isnonzero[1:T] .* meas_conc[1:T];
    predicted_concentration_norm = predicted_concentration .* flow[1:T];
    if (h>0) {
      predicted_concentration_forecast = isnonzero[(T+1):(T+h)] .* meas_conc[(T+1):(T+h)];
      predicted_concentration_forecast_norm = predicted_concentration_forecast .* flow[(T+1):(T+h)];
    }
  }
}
