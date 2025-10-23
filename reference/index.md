# Package index

## EpiSewer

- [`sewer_data()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_data.md)
  : Specify observation data
- [`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md)
  : Specify modeling assumptions
- [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  : Estimate epidemiological parameters from wastewater measurements

## Measurements module

Modeling functions for the measurement process

- [`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md)
  : Model the measurement process
- [`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md)
  : Observe concentration measurements
- [`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md)
  : Estimate measurement noise
- [`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md)
  : Estimate measurement noise with constant variance
- [`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
  : Estimate measurement noise for digital PCR data
- [`LOD_assume()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_assume.md)
  : Assume a limit of detection
- [`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
  : Estimate a limit of detection model for digital PCR data
- [`LOD_none()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_none.md)
  : Do not model a limit of detection

## Sampling module

Modeling functions for the sampling process

- [`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md)
  : Model the sampling process
- [`outliers_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/outliers_estimate.md)
  : Model outliers via an extreme value distribution
- [`outliers_none()`](https://adrian-lison.github.io/EpiSewer/reference/outliers_none.md)
  : Do not model outliers
- [`sample_effects_estimate_matrix()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_matrix.md)
  : Estimate sample effects using a design matrix
- [`sample_effects_estimate_weekday()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_weekday.md)
  : Estimate weekday sample effects
- [`sample_effects_none()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_none.md)
  : Do not model sample effects

## Sewage module

Modeling functions for the sewage process

- [`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md)
  : Model the sewage process
- [`flows_assume()`](https://adrian-lison.github.io/EpiSewer/reference/flows_assume.md)
  : Assume a constant wastewater flow
- [`flows_observe()`](https://adrian-lison.github.io/EpiSewer/reference/flows_observe.md)
  : Observe wastewater flows
- [`residence_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/residence_dist_assume.md)
  : Assume a sewer residence time distribution

## Shedding module

Modeling functions for the shedding process

- [`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md)
  : Model the shedding process
- [`shedding_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/shedding_dist_assume.md)
  : Assume a shedding load distribution
- [`shedding_dist_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/shedding_dist_estimate.md)
  : Estimate an uncertain shedding load distribution
- [`incubation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/incubation_dist_assume.md)
  : Assume an incubation period distribution
- [`load_per_case_assume()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_assume.md)
  : Assume the average load per case
- [`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md)
  : Calibrate the average load per case using case count data
- [`load_variation_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/load_variation_estimate.md)
  : Estimate individual-level load variation
- [`load_variation_none()`](https://adrian-lison.github.io/EpiSewer/reference/load_variation_none.md)
  : Do not model individual-level load variation

## Infections module

Modeling functions for the infection process

- [`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md)
  : Model the infection process
- [`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md)
  : Estimate Rt using Gaussian processes
- [`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
  : Estimate Rt via smoothing splines
- [`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md)
  : Estimate Rt via a random walk
- [`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md)
  : Estimate Rt via exponential smoothing
- [`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md)
  : Estimate Rt via piecewise constant changepoint model
- [`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md)
  : Estimate Rt with a smooth derivative
- [`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md)
  : Estimate Rt via a changepoint spline model
- [`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md)
  **\[deprecated\]** : Estimate Rt using an approximation of the
  generative renewal model
- [`generation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/generation_dist_assume.md)
  : Assume a generation time distribution
- [`seeding_estimate_constant()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_constant.md)
  : Estimate constant seeding infections
- [`seeding_estimate_growth()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_growth.md)
  : Estimate seeding infections with a time-varying growth rate
- [`seeding_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_rw.md)
  : Estimate seeding infections using a random walk model
- [`infection_noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/infection_noise_estimate.md)
  : Estimate infection noise
- [`infection_noise_none()`](https://adrian-lison.github.io/EpiSewer/reference/infection_noise_none.md)
  : Do not model infection noise

## Forecast module

Modeling functions for specifying forecasts

- [`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md)
  : Forecasting module
- [`horizon_assume()`](https://adrian-lison.github.io/EpiSewer/reference/horizon_assume.md)
  : Specify the forecast horizon
- [`horizon_none()`](https://adrian-lison.github.io/EpiSewer/reference/horizon_none.md)
  : Do not produce forecasts
- [`damping_assume()`](https://adrian-lison.github.io/EpiSewer/reference/damping_assume.md)
  : Dampen forecasts
- [`damping_none()`](https://adrian-lison.github.io/EpiSewer/reference/damping_none.md)
  : Do not dampen forecasts

## Model fitting

Functions to specify and run the sampling of EpiSewer models

- [`set_fit_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_fit_opts.md)
  : Configure the model fitting
- [`sampler_stan_mcmc()`](https://adrian-lison.github.io/EpiSewer/reference/sampler_stan_mcmc.md)
  : Use the stan MCMC sampler
- [`sampler_stan_pathfinder()`](https://adrian-lison.github.io/EpiSewer/reference/sampler_stan_pathfinder.md)
  : Use stan's pathfinder variational inference algorithm
- [`model_stan_opts()`](https://adrian-lison.github.io/EpiSewer/reference/model_stan_opts.md)
  : Specify details of the stan model
- [`set_results_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_results_opts.md)
  : Configure results returned after model fitting
- [`sewer_compile()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_compile.md)
  : Compile EpiSewer models
- [`run()`](https://adrian-lison.github.io/EpiSewer/reference/run.md) :
  Fit an EpiSewer model.

## Plotting

Functions for plotting

- [`plot_R()`](https://adrian-lison.github.io/EpiSewer/reference/plot_R.md)
  : Plot the effective reproduction number
- [`plot_growth_report()`](https://adrian-lison.github.io/EpiSewer/reference/plot_growth_report.md)
  : Plot a growth report
- [`plot_growth_rate()`](https://adrian-lison.github.io/EpiSewer/reference/plot_growth_rate.md)
  : Plot the epidemic growth rate
- [`plot_doubling_time()`](https://adrian-lison.github.io/EpiSewer/reference/plot_doubling_time.md)
  : Plot the epidemic doubling time
- [`plot_infections()`](https://adrian-lison.github.io/EpiSewer/reference/plot_infections.md)
  : Plot infections
- [`plot_load()`](https://adrian-lison.github.io/EpiSewer/reference/plot_load.md)
  : Plot the estimated load
- [`plot_concentration()`](https://adrian-lison.github.io/EpiSewer/reference/plot_concentration.md)
  : Plot predicted concentration
- [`plot_prior_posterior()`](https://adrian-lison.github.io/EpiSewer/reference/plot_prior_posterior.md)
  : Visually compare prior and posterior of a model parameter
- [`plot_LOD()`](https://adrian-lison.github.io/EpiSewer/reference/plot_LOD.md)
  : Plot limit of detection
- [`plot_sample_effects()`](https://adrian-lison.github.io/EpiSewer/reference/plot_sample_effects.md)
  : Plot estimated sample effects

## Specification

Helper functions for specifying distributions and other assumptions

- [`get_discrete_gamma()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_gamma.md)
  : Get PMF of a discretized Gamma distribution
- [`get_discrete_gamma_shifted()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_gamma_shifted.md)
  : Get PMF of a discretized shifted Gamma distribution
- [`get_discrete_lognormal()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_lognormal.md)
  : Get PMF of a discretized lognormal distribution
- [`suggest_load_per_case()`](https://adrian-lison.github.io/EpiSewer/reference/suggest_load_per_case.md)
  : Suggest load per case assumption using wastewater data and case
  numbers
- [`qexpgamma()`](https://adrian-lison.github.io/EpiSewer/reference/qexpgamma.md)
  : Exponential-Gamma distribution

## Preprocessing

Functions for preprocessing data

- [`mark_outlier_spikes_median()`](https://adrian-lison.github.io/EpiSewer/reference/mark_outlier_spikes_median.md)
  **\[experimental\]** : Mark outlier spikes in a measurement time
  series

## Utilities

Utility functions

- [`component_functions()`](https://adrian-lison.github.io/EpiSewer/reference/component_functions.md)
  : Get a list of modeling functions for a component
- [`modeldata_init()`](https://adrian-lison.github.io/EpiSewer/reference/modeldata_init.md)
  : Construct an unspecified EpiSewer model
- [`get_checksums()`](https://adrian-lison.github.io/EpiSewer/reference/get_checksums.md)
  : Get checksums that uniquely identify an EpiSewer job.
- [`summarize_fit()`](https://adrian-lison.github.io/EpiSewer/reference/summarize_fit.md)
  : Summarize parameters of interest
