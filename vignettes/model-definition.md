# EpiSewer model definition
Adrian Lison\
ETH Zurich, Department of Biosystems Science and Engineering, Zurich, Switzerland

This document details the generative model used by EpiSewer to estimate reproduction numbers from wastewater concentration measurements. The model is structured into 5 modules, namely [infections](#infections-module), [shedding](#shedding-module), [sewage](#sewage-module), [sampling](#sampling-module), and [measurements](#measurements-module).

Note that below, the full model with all options is specified, while `EpiSewer` excludes certain parts of the model by default (e.g. no multi-day composite samples) - these have to be explicitly selected by the user to be included.
## Infections module
#### Reproduction number ($R_t$)
##### Link function
The model employs a link function to ensure that $R_t>0$ regardless of the smoothing approach used. An obvious choice for this is the log link, which however has the disadvantage that it implies an asymmetric prior on changes in $R_t$ (for example, a change from $R_t=2$ to $R_t=1$ would be much more likely than a change from $R_t=1$ to $R_t=2$). `EpiSewer` therefore provides alternative link functions which are configured to behave approximately like the identity function around $R_t=1$ but ensure that $R_t>0$.
###### Option 1: Inverse softplus
$$\text{link}^{-1}(x) = \text{softplus}_k(x) = \frac{\log(1+e^{kx})}{k}$$
with sharpness parameter $k$.

By default, we choose $k=4$, which makes the softplus function behave effectively like $f(x)=x$ for practically relevant values of $R_t$ while ensuring that $R_t>0$.
###### Option 2: Scaled logit 
$$\text{link}^{-1}(x) = \text{logistic}_{c, a, k}(x) = \frac{c}{1 + a \, e^{-k \, x}}$$
with hyperparameters $c$, $a$, and $k$. 

Here, $c$ defines the maximum possible $R_t$ (default $R_\max=6$). This link function is inspired from the epidemia package[^epidemia], but has two further hyperparameters $a$ and $k$, which are chosen such that the logistic function equals 1 at $x=1$ and also has a derivative of 1. Thus, the link function behaves roughly like the identity $f(x)=x$ around $R_t=1$, which makes the definition of adequate smoothing priors for $R_t$ easier.
##### Smoothing
$R_t$ is modeled via a time series or smoothing spline prior which ensure smoothness through temporal autocorrelation of the $R_t$ estimates.
###### Option 1: Random walk
$$\text{link}(R_t)\,|\,R_{t-1} \sim N(\text{link}(R_{t-1}),{\sigma_{\text{link}(R)}}^2)$$
with a normal prior on the intercept $R_1$ and a zero-truncated normal prior on the standard deviation $\sigma_{\text{link}(R)}$ of the random walk.
###### Option 2: Exponential smoothing
This option uses an **innovations state space model** with additive errors to implement an exponential smoothing prior (Holt's linear trend method with dampening)[^hyndman]. For this, $R_t$ is modeled recursively as
$$\begin{align*}
    \text{link}(R_t) &= l_{t-1} + \phi' \, b_{t-1} + \epsilon_t \\
    l_t &= l_{t-1} + \phi' \, b_{t-1} + \alpha' \, \epsilon_t \\
    b_t &= \phi' \, b_{t-1} + \beta' \, \epsilon_t \\
    \epsilon_t &\sim N(0,{\sigma_{\text{link}(R)}}^2),
\end{align*}$$
where $l_t$ is the level, $b_t$ the trend, and $\epsilon_t$ the additive error at time $t$, respectively. 

Moreover, $\alpha'$ is a level smoothing parameter, $\beta'$ a trend smoothing parameter, and $\phi'$ a dampening parameter. We use normal priors for the intercepts of the level $l_1$ and trend $b_1$, and a zero-truncated normal prior for the variance ${\sigma_{\text{link}(R)}}^2$ of the innovations, allowing these to be estimated from the data. We use Beta priors for the smoothing parameters $\alpha'$ and $\beta'$.

In theory, the dampening parameter $\phi'$ can also be estimated from the data, but usually has to be fixed due to identifiability issues. The default is $\phi'=0.9$, i.e. the strength of the trend is assumed to halve approximately every week of predicting into the future.
###### Option 3: Penalized B-splines
$$\text{link}(R_t) = \boldsymbol{B}_{[t,\,]} \,  a $$
where $\boldsymbol{B}_{[t,\,]}$ is a (by default cubic) B-spline basis evaluated at date index $t$, and $a = (a_1, a_2, \dots, a_\text{df})$ is a vector of spline coefficients, with $\text{df} = \text{number of knots} + \text{spline degree}$ being the degrees of freedom.

We here use *penalized* B-splines, which we implemented by regularizing the spline coefficients via a random walk[^Kharratzadeh], i.e.
$$a_{i}\,|\,a_{i-1} \sim N(a_{i-1}, \sigma_{a}^2 \, \delta_i)$$
where the random walk variance $\sigma_{a}^2$ is scaled by $\delta_i$, the distance between knot $i$ and knot $i-1$. This allows to provide a (zero-truncated normal) prior on $\sigma_{a}$ which is independent of the distance between knots. By enforcing smoothness of the spline coefficients, the risk of overfitting the splines is greatly reduced. For the intercept $a_1$, a normal prior reflecting the expectation of $R_t$ at the beginning of the time series is used.

Knots are placed uniformly according to a user-specified knot distance, however the knot distance is heuristically reduced for dates close to the present (which are only partially informed by data) to avoid erroneous extrapolation.
#### Infections
EpiSewer uses a stochastic renewal model as described in earlier work[^semimechanistic][^generative_nowcasting] to model both expected and realized infections.
##### Expected infections ($\iota_t$)
###### Seeding
We require a seeding phase at the start of the modeled time period, when the renewal model cannot be applied yet. Here we use a random walk on the log scale to model the expected number of initial infections.
$$\log(\iota_t)|\iota_{t-1} \sim N(\log(\iota_{t-1}),{\sigma_{\log(\iota)}}^2) \;|\; t \leq G$$
with a normal prior on intercept $\log(\iota_1)$ and a zero-truncated normal prior on the standard deviation $\sigma_{\log(\iota)}$ of the random walk.
###### Renewal model

$$ \iota_t = E[I_t] = R_t \sum_{s=1}^G \tau^\text{gen}_s \, I_{t-s} \;|\; t > G $$
where $I_{t-s}$ is the number of realized infections that occurred $s$ days before $t$, and $\tau^\text{gen} = (\tau^\text{gen}_1, \tau^\text{gen}_2, \dots, \tau^\text{gen}_G)$ the discrete generation time distribution with maximum generation time $G$.
##### Realized infections ($I_t$)
###### Option 1: Poisson
$$I_t|\iota_t \sim \text{Poisson}(\iota_t)$$
with expected infections $\iota_t$.

In the stan implementation, this is approximated as $I_t|\iota_t \sim N\left(\iota_t, \iota_t\right)$
###### Option 2: Negative binomial
$$I_t|\iota_t \sim \text{Negative binomial}(\mu = \iota_t, \sigma^2 = \iota_t + \frac{\iota_t^2}{\phi})$$
with expected infections $\iota_t$ and inverse overdispersion parameter $\phi$ . Following  [^prior_recommendations], we place a zero-truncated normal prior not on $\phi$ but rather on $\xi = \frac{1}{\sqrt{\phi}}$. 

In the stan implementation, this is approximated as $I_t|\iota_t \sim N\left(\iota_t, \iota_t + \frac{\iota_t^2}{\phi}\right)$
## Shedding module
##### Expected symptom onsets ($\lambda_t$)
$$ \lambda_t = \sum_{s=0}^L I_{t-s}\, \tau^\text{inc}_s$$
where $\tau^\text{inc} = (\tau^\text{inc}_0, \tau^\text{inc}_1, \dots, \tau^\text{inc}_L)$ is the incubation period distribution with maximum incubation period $L$.
##### Total load shed in catchment ($\omega_t$)
###### Option 1: Without individual-level variation

$$ \omega_t = \sum_{s=0}^S \zeta_{t-s}\, \tau^\text{shed}_s$$
where  $\zeta_{t} = \mu^\text{load}\, \lambda_{t}$  is the total load shed on day $t$, with an average total shedding load per person of $\mu^\text{load}$, and $\tau^\text{shed} = (\tau^\text{shed}_0, \tau^\text{shed}_1, \dots, \tau^\text{shed}_S)$ is the distribution of shedding over time, with shedding up to $S$ days after symptom onset. 
###### Option 2: With individual-level variation

$$ \omega_t = \sum_{s=0}^S \zeta_{t-s}\, \tau^\text{shed}_s$$
where $\tau^\text{shed} = (\tau^\text{shed}_0, \tau^\text{shed}_1, \dots, \tau^\text{shed}_S)$ is again the distribution of shedding over time, with shedding up to $S$ days after symptom onset, but $\zeta_t$ is now defined as the sum of i.i.d. gamma distributed individual shedding loads with mean $\mu^\text{load}$ and coefficient of variation  $\nu_\zeta$:
$$\zeta_t \sim \sum_{i=1}^{\Lambda_t} \, \text{Gamma}(\text{mean} = \mu^\text{load}, \text{cv} = \nu_\zeta) = \text{Gamma}(\alpha = \frac{\Lambda_t}{\nu_\zeta^2}, \beta = \frac{1}{\mu^\text{load}\, \nu_\zeta^2}) $$
where $\Lambda_t$ is the realized number of symptom onsets. We here approximate $\Lambda_t$ with its expectation $\lambda_t = E[\Lambda_t]$, i.e.
$$\zeta_t \sim \text{Gamma}(\alpha = \frac{\lambda_t}{\nu_\zeta^2}, \beta = \frac{1}{\mu^\text{load}\, \nu_\zeta^2}),$$
and place a zero-truncated normal prior on $\nu_\zeta$.
## Sewage module
#### Load at sampling site ($\pi_t$)

$\pi_t =  \sum_{d=0}^D \omega_{t-d}\, \tau^\text{res}_d$

where $\tau^\text{res} = (\tau^\text{res}_0, \tau^\text{res}_1, \dots, \tau^\text{res}_D)$ is the sewer residence time distribution of pathogen particles.
#### Concentration at sampling site ($\chi_t$)

$$\chi_t = \frac{\pi_t}{\text{flow}_t}$$
where $\text{flow}_t$ is a measure of the total wastewater volume upstream from the sampling site on day $t$.
## Sampling module
#### Concentration in single-day samples ($\kappa_t$), with sample effects
$$ \kappa_t = \chi_t + \boldsymbol{X} \, \eta $$
where $\boldsymbol{X}$ is a $T \times K$ design matrix of $K$ different covariates describing the characteristics of samples on days 1 to $T$, and $\eta = (\eta_1, \eta_2, \dots, \eta_K)$ is a vector of coefficients which estimate the association of the sample covariates with the concentration. We place independent normal priors on $\eta_1, \eta_2, \dots, \eta_K$.

Note that if on some days no samples were taken, the corresponding row of the design matrix $\boldsymbol{X}$ can be filled with arbitrary values, since they will not contribute to the likelihood.
#### Concentration in multi-day composite samples ($\rho_t$)
$$\rho_t = \frac{1}{w} \sum_{i=0}^{w-1} \kappa_{t-i}$$
where $w$ is the composite window size, i.e. the number of consecutive single-day samples that are equivolumetrically mixed into a composite sample. Here, $t$ is the day of the latest sample.
## Measurements module
#### Analysed concentration before replication stage ($\Psi_t$)
We model $\Psi_t$ as Log-Normal distributed with unit-scale mean $\rho_t$ (sample concentration) and unit-scale coefficient of variation $\nu_\psi$, i.e.
$$ \Psi_t \sim \text{Log-Normal}\left(\mu_\Psi,\, \sigma_\Psi^2\right)$$
where $\mu_\Psi = \log(\rho_t) - \frac{1}{2}\sigma_\Psi^2$ is the location and $\sigma_\Psi^2 = \log(1 + \nu_\Psi^2)$ is the scale of the Log-Normal distribution. Here, $\nu_\psi$ measures the unexplained variation in concentrations before the replication stage (for example, if the replication stage is PCR quantification, then $\nu_\psi$ measures all variation before PCR). We place a zero-truncated normal prior on $\nu_\psi$.
#### Limit of detection (LOD)
We implement a hurdle model for the limit of detection, where the probability for a zero measurement is described by a logistic regression model with midpoint at the LOD, i.e. $$ P(\Upsilon_{t,i} = 0 | \Psi_t = \psi_t) = \frac{1}{1 + e^{-\frac{\text{LOD} - \psi_t}{\text{LOD}}\, k}}$$where $\Upsilon_{t,i}$ is the $i^{\text{th}}$ replicate concentration measurement from the (composite) sample on day $t$. Moreover, $k$ is a sharpness parameter describing the steepness of the logistic function (the steeper the function, the sharper the modeled limit of detection). For non-zero measurements, we conversely have
$$P(\Upsilon_{t,i} > 0 | \Psi_t = \psi_t) = 1 - P(\Upsilon_t = 0 | \Psi_t = \psi_t)$$
#### Measured concentration ($\Upsilon_t$)
Non-zero concentration measurements are modeled as log-normal distributed, i.e.
$$ \Upsilon_{t,i} \,|\, \Upsilon_{t,i} > 0 \, \sim \text{Log-Normal}\left(\mu_{\Upsilon_{t}},\, \sigma_\Upsilon^2\right)$$
where $\Upsilon_{t,i}$ is the $i^{\text{th}}$ replicate concentration measurement from the (composite) sample on day $t$. Here, $\mu_{\Upsilon_{t}} = \log(\psi_t) - \frac{1}{2}\sigma_{\Upsilon}^2$ is the location and $\sigma_\Upsilon^2 = \log(1 + \nu_\Upsilon^2)$ the scale of the Log-Normal distribution, with $\nu_\Upsilon$ being the coefficient of variation for the replication stage. This means that if, for example, the replication stage is PCR quantification, then $\nu_\Upsilon$ measures the noise from the PCR. We place a zero-truncated normal prior on $\nu_\Upsilon$.

Note that if no replicates are modeled, then $\nu_\Upsilon$ measures the total variation and $\nu_\psi$ is not estimated.
## Estimation
EpiSewer is directly fitted to observed concentration measurements $y_{t,i}$ (where $t<T$ is the index of the sample date and $i$ the index of the replicate out of $n_t$ replicates in total for date index $t$).

Using Markov Chain Monte Carlo (MCMC), we draw samples from the joint posterior distribution of $R_t$ and the other parameters $$P\left(R_{\{t | t \leq T\}}, \dots \;\middle|\; y_{\{t | t \leq T\},\{i | i \leq n_t\}}\right) \propto P\left(y_{\{t | t \leq T\},\{i | i \leq n_t\}} \;\middle|\; R_{\{t | t \leq T\}}, \dots\right) \, P\left(\dots \right),$$where $P\left(y_{\{t | t \leq T\},\{i | i \leq n_t\}} \;\middle|\; R_{\{t | t \leq T\}}, \dots\right)$ is the likelihood as defined by the 5 modules of our generative model (infection, shedding, sewage, sampling, and measurement), and $P\left(\dots \right)$ represents the priors on the non-hierarchical parameters of the model (e.g. intercepts and variances of random walks).

Note that many variables in the generative model (including $R_t$!) are in fact only transformations of other variables and therefore do not count as model parameters in a strict sense.
## Indexing of time
For ease of notation, we have here used a unified time index $t \in \{1, 2, \dots, T\}$ for all variables. In reality, to model $\Upsilon_{t,i}$ at $t=1$, we need to model $R_t$ and many other variables also at earlier time points $t<1$ due to the various delays between infections and observations. In the model implementation, each variable therefore has its own time index, and indices are mapped to each other accordingly.
## References

[^epidemia]: Scott, J. A. _et al._ epidemia: Modeling of epidemics using hierarchical bayesian models. (2020). Available from: https://imperialcollegelondon.github.io/epidemia/index.html
[^hyndman]: Hyndman, R. J. & Athanasopoulos, G. _Forecasting: principles and practice_. (OTexts, 2018).
[^Kharratzadeh]: Kharratzadeh, M. & Stan development team. Splines in Stan. Stan documentation case study. (2017). Available from: https://mc-stan.org.
[^prior_recommendations]: Stan development team. Prior Choice Recommendations. (2020). Available from: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations.
[^generative_nowcasting]: Lison, A., Abbott, S., Huisman, J. & Stadler, T. Generative Bayesian modeling to nowcast the effective reproduction number from line list data with missing symptom onset dates. Preprint at [https://doi.org/10.48550/arXiv.2308.13262](https://doi.org/10.48550/arXiv.2308.13262) (2023).
[^semimechanistic]: Bhatt, S. _et al._ Semi-Mechanistic Bayesian Modeling of COVID-19 with Renewal Processes. _arXiv:2012.00394_ (2020).