/**
  * Compute coefficient of variation (CV) in ddPCR quantification as a function
  * of the concentration
  *
  * @param concentration vector with concentrations (typically gc/mlWW)
  *
  * @param a intercept of the CV (all other independent noise factors)
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
  return(a + sqrt(expm1(conc_scaled)/b) ./ conc_scaled);
}
