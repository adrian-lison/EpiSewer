/**
  * Compute coefficient of variation (CV) in ddPCR quantification as a function
  * of the concentration
  *
  * @param concentration vector with concentrations (typically gc/mlWW)
  *
  * @param a pre-PCR CV (all other noise factors)
  *
  * @param b number of droplets
  *
  * @param c scaling factor, should be $c = \text{droplet volume} / (\frac{\text{reaction volume}}{\text{template volume}}\times \text{dilution} \times \frac{\text{elution volume}}{\text{volume mL}} \times 1000)$
  *
  * @return A vector with the corresponding coefficients of variation
  */
vector cv_ddPCR(vector concentration, real a, real b, real c) {
  int n = num_elements(concentration);
  vector[n] conc_scaled = concentration * c;
  vector[n] pre_PCR_factor = sqrt(1 + a^2 * b * conc_scaled);
  return(pre_PCR_factor .* sqrt(expm1(conc_scaled)/b) ./ conc_scaled);
}
