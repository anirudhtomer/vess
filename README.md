# VESS: An R package for vaccine efficacy sample size calculations in a multiple vaccine and multiple pathogen variant scenario.

## Motivation
An important problem in vaccine efficacy (VE) studies is sample size calculation during study planning. Vaccine manufacturers may face two types of sample size challenges. First, to have enough subjects to show a non-zero VE or to show that VE is above a certain threshold (superiority trial). In this setting one can obtain an appropriate sample size utilizing the hypothesis testing framework. A second type of challenge is having enough subjects to obtain precise estimates of VE to further assist in decisions pertaining to the various phases of the development of a vaccine. In statistical terms, precision refers to confidence limits of the VE, and a narrower confidence interval can be obtained with a larger sample size. While the actual confidence limits depend upon the actual data, using Monte Carlo simulations confidence intervals can also be simulated and a range for both the upper and lower confidence limit can be obtained. Currently available tools such as the **R** package **epiR** (<https://cran.r-project.org/web/packages/epiR/index.html>) and general purpose calculators such as PASS (<https://www.ncss.com/software/pass/>) only handle the single variant and single vaccine scenario. Another drawback of the current tools is that for sample size calculations aiming at obtaining a certain precision (Kelley et al., 2003) for VE, the calculators provide precision assuming that the actual data will be equal to the expected cell counts of a 2 x 2 cross table of infection and vaccination status. Thus, they ignore the randomness of the data generating process. Besides, current calculators also do not support multiple variants and multiple vaccines.

## Features
This package has 3 R functions for finding the minimum detectable vaccine efficacy given a sample size, power, and Type-I error. 
* `get_cohort_mindet_VE_irr()`
* `get_cohort_mindet_VE_cir()`
* `get_cohort_mindet_VE_or()`

Then there are 3 R functions for finding the expected precision for a given vaccine efficacy, sample size, and Type-I error.
* `get_cohort_expectedCI_VE_irr()`
* `get_cohort_expectedCI_VE_cir()`
* `get_cohort_expectedCI_VE_or()`

Lastly, there are 2 R functions for 
* `get_allornone_ve_estimands()`
* `get_leaky_ve_estimands()`

## References
Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or obtaining precision: Delineating methods of sample-size planning. Evaluation & the health professions, 26(3), 258-287.

