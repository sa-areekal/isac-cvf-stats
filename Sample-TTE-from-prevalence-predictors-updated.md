Sample size calculations for time-to-events analysis prevalence.
================

``` r
calc_ss_tte_crit1 <- function(p, R2_CS_adj, S_VH = 0.9) {
  
  
  # Calculate the sample size n
  n <- p / ((S_VH - 1) * log(1 - (R2_CS_adj / S_VH)))
  
  
  # Return the calculated sample size
  return(n)
}
```

# STRATIFY

### Stratify calculations

Calculate likelihood of the null model $$
\ln L_{\text{null}} = E \ln(\lambda) + \lambda T = E \ln\left(\frac{E}{T}\right) + E
$$

``` r
# event rate of 24.6 per 100000 person year, median follow-up of 7 years
# (log_L_null <- N_VD*log(N_VD/total_person_year) + N_VD)
(stratify_lnLnull <- (24.6 * log(24.6/100000)) + 24.6 )
```

    ## [1] -179.8304

``` r
(stratify_lnLnull_logreg <- (24.6 * log(24.6/100000/7)) + (100000/7 - 24.6)*log(1-(24.6/100000/7)))
```

    ## [1] -252.801

max R^2\_{app}

``` r
N_personyear <- 100000
N_tot_stratify <- 100000/7 # total person-years/follow-up time
(stratify_max_R2_cs_app <- 1 - exp(2*stratify_lnLnull/N_tot_stratify))
```

    ## [1] 0.02486198

``` r
# (stratify_max_R2_cs_app <- 1 - exp(2*stratify_lnLnull_logreg/N_tot_stratify))
#(stratify_max_R2_cs_app <- 1 - exp(2*stratify_lnLnull/N_personyear))
```

R^2_CS_adj

$$
R^2_{CS\_app} = R^2_{Nagelkerke\_app} \left(1 - \exp\left(\frac{2 \ln L_{\text{null}}}{n}\right)\right) 
$$

``` r
(stratify_R_2_app_new <- 0.15*(1-exp(2*stratify_lnLnull/N_tot_stratify)))
```

    ## [1] 0.003729297

For Nagelkerke’s R-squared statistic to be 0.15, the $R^2_{app}$ needs
to be

``` r
(stratify_R_2_app <-  0.15 * stratify_max_R2_cs_app )
```

    ## [1] 0.003729297

From Eq. 8, the estimated $R^2_{adj}$ is

``` r
(stratify_R_2_adj_cal <- 0.9 * stratify_R_2_app)
```

    ## [1] 0.003356367

S_VH based on crit 2

``` r
calc_S_VH <- function(R2_CS_adj, delta, max_R2_CS_app) {
  # Calculate S_VH using the given formula
  S_VH <- R2_CS_adj / (R2_CS_adj + delta * max_R2_CS_app)
  
  # Return the calculated S_VH
  return(S_VH)
}

calc_S_VH(R2_CS_adj = 0.003, delta = 0.05, max_R2_CS_app = stratify_max_R2_cs_app)
```

    ## [1] 0.7070304

``` r
(calc_ss_tte_crit1(p = 40, R2_CS_adj = stratify_R_2_adj_cal, S_VH = 0.90))
```

    ## [1] 107058.7

``` r
(calc_ss_tte_crit1(p = 40, R2_CS_adj = 0.004, S_VH = 0.90))
```

    ## [1] 89799.85

# CVF Phenotype

## Prevalence of outcomes under study.

Vacular dementia (VD) incidence from individuals with hypertension
(Adesuyan et al. 2023)

N = 85,245 and

VD = 19,581

### From Chung et al. (2023):

Total number of eligible individuals: 4,309,481; Total number of
incident dementia patients: 145,319 (3.4% of the population); Percentage
of dementia patients with vascular dementia: 46%; Median follow-up time:
9.5 years

The number of vascular dementia cases can be calculated as:

``` r
N_tot <- 4309481 # total no of subjects
N_dem <- 145319 # number of dementia cases
t_foll <- 9.5 # years
(N_VD <- 0.46 * N_dem) # Vascular dementia cases
```

    ## [1] 66846.74

The incidence rate is typically expressed as the number of cases per
person-year of follow-up

``` r
total_person_year <- N_tot * t_foll
print(paste("Total person years is",total_person_year))
```

    ## [1] "Total person years is 40940069.5"

``` r
VD_inc <- N_VD/total_person_year
print(paste("Vascular dementia incidence is",VD_inc*100000, "cases per 100000 persons"))
```

    ## [1] "Vascular dementia incidence is 163.279498096602 cases per 100000 persons"

Assuming an event rate of 163 per 100,000 years, a median follow up of
9.5 years, a conservative estimate of Nagelkerke’s R-squared statistic
of 0.15 and a maximum number of 40 parameters per model, the number of
sample size required is:

### Calculating max($R^2_{CS_app}$)

For example, assuming (for simplicity) an exponential survival model was
fitted

``` r
(log_L_null <- N_VD*log(N_VD/total_person_year) + N_VD)
```

    ## [1] -362139.7

$\max{R^2_{app}}$ will be

``` r
(max_R2_cs_app <- 1 - exp(2*log_L_null/N_tot))
```

    ## [1] 0.1547024

For Nagelkerke’s R-squared statistic to be 0.15, the $R^2_{app}$

``` r
(R_2_app <-  0.15 * max_R2_cs_app )
```

    ## [1] 0.02320536

From Eq. 8, the estimated $R^2_{adl}$ is

``` r
(R_2_adj_cal <- 0.9 * R_2_app)
```

    ## [1] 0.02088482

``` r
p <- 35 # number of predictor parameters


calc_ss_tte_crit1 <- function(p, R2_CS_adj, S_VH = 0.9) {
  
  
  # Calculate the sample size n
  n <- p / ((S_VH - 1) * log(1 - (R2_CS_adj / S_VH)))
  
  
  # Return the calculated sample size
  return(n)
}
```

New sample size calculation for Vascular dementia

``` r
calc_ss_tte_crit1(p = 30, R2_CS_adj = R_2_adj_cal, S_VH = 0.90)
```

    ## [1] 12777.46

S_VH based on crit 2

``` r
calc_S_VH <- function(R2_CS_adj, delta, max_R2_CS_app) {
  # Calculate S_VH using the given formula
  S_VH <- R2_CS_adj / (R2_CS_adj + delta * max_R2_CS_app)
  
  # Return the calculated S_VH
  return(S_VH)
}

 calc_S_VH(R2_CS_adj = R_2_adj_cal, delta = 0.05, max_R2_CS_app = stratify_max_R2_cs_app)
```

    ## [1] 0.9438222

``` r
(sample_VaD <- calc_ss_tte_crit1(p = 26, R2_CS_adj = R_2_adj_cal, S_VH = 0.94))
```

    ## [1] 19286.32

## crit 3

``` r
1-exp(-0.001632795*9.5)
```

    ## [1] 0.01539187

``` r
# apprx 29671.26 as 30,000
(1-exp(-(0.00163+1.96*sqrt(0.00163/30000*9.5))*9.5))-
(1-exp(-(0.00163-1.96*sqrt(0.00163/30000*9.5))*9.5))
```

    ## [1] 0.02634466

less that the accepted value of delta 0.05

Samples required

``` r
# 30000 participants
print(paste("The number of patients required is", round(sample_VaD)))
```

    ## [1] "The number of patients required is 19286"

``` r
# person-years of follow-up
print(paste("The person years of follow up", round(sample_VaD*9.5)))
```

    ## [1] "The person years of follow up 183220"

``` r
print(paste("The number of outcome events expected is", round(sample_VaD*9.5*VD_inc)))
```

    ## [1] "The number of outcome events expected is 299"

``` r
print(paste("EPP expected is", round(sample_VaD*9.5*VD_inc/40)))
```

    ## [1] "EPP expected is 7"

For 100 EPP, sample size required is

``` r
100*40/(9.5*VD_inc)
```

    ## [1] 257872.3

### References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-adesuyan2023" class="csl-entry">

Adesuyan, Matthew, Yogini H. Jani, Dana Alsugeir, Robert Howard, Ian C.
K. Wong, Li Wei, and Ruth Brauer. 2023. “Trends in the Incidence of
Dementia in People with Hypertension in the UK 2000 to 2021.”
*Alzheimer’s & Dementia : Diagnosis, Assessment & Disease Monitoring* 15
(3): e12466. <https://doi.org/10.1002/dad2.12466>.

</div>

<div id="ref-chung2023" class="csl-entry">

Chung, Sheng-Chia, Rui Providencia, Reecha Sofat, Mar Pujades-Rodriguez,
Ana Torralbo, Ghazaleh Fatemifar, Natalie K. Fitzpatrick, et al. 2023.
“Incidence, Morbidity, Mortality and Disparities in Dementia: A
Population Linked Electronic Health Records Study of 4.3 Million
Individuals.” *Alzheimer’s & Dementia* 19 (1): 123–35.
<https://doi.org/10.1002/alz.12635>.

</div>

</div>
