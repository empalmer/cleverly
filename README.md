
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cleverly

<!-- badges: start -->
<!-- badges: end -->

Cleverly stands for (Cl)ustering with (E)xternal (V)ariabl(e)s in (R)
with (L)ongitudinal (Y)s.

## Installation

You can install the development version of cleverly from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("empalmer/cleverly")
library(cleverly)
```

``` r
set.seed(127)
  sim <- sim_Z_longitudinal(n = 20,
                            range_start = 5000,
                            range_end = 20000,
                            nknots = 3,
                            K = 12,
                            order = 3,
                            user_var = 1000,
                            cor_str = "IND",
                            al = 0.4)

  Y <- dplyr::select(sim, -c(
                        "total_n",
                        "Capture.Number",
                        "Z"))
  Z <- sim$Z

  start <- Sys.time()
  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  lp = 0,
                  time = time,
                  gammas = c(1, 1), # controls smoothness
                  tau = 0.1, # Controls cuttoff for highest shrinkage
                  theta = 3000, # for lambda, but also for d
                  psi = 2900, # controls clustering
                  C = 100,
                  max_admm_iter = 50,
                  max_outer_iter = 5,
                  max_2_iter = 100,
                  epsilon_r = .001,
                  epsilon_d = .05,
                  epsilon_b = .01,
                  epsilon_2 = .001,
                  cor_str = "IND")
  end <- Sys.time()
  (duration <- end - start)
  res$clusters$no
```

<!--
This is a basic example which shows you how to solve a common problem:
&#10;
``` r
library(cleverly)
## basic example code
```
&#10;What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:
&#10;
``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```
&#10;You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.
&#10;You can also embed plots, for example:
&#10;<img src="man/figures/README-pressure-1.png" width="100%" />
&#10;In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
-->
