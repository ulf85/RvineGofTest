# R Code for "Generalized Information Matrix Tests for Copulas"
In this repository I publish my **R code** for the paper **"Generalized Information Matrix Tests for Copulas"** written with my colleagues Artem Prokhorov and Yajing Zhu (not yet published but accepted in Econometric Reviews).

A copy is available on my [researchgate](https://www.researchgate.net/publication/299641789_Generalized_Information_Matrix_Tests_for_Copulas) page.

In this paper we propose a family of goodness-of-fit tests for copulas. The tests use generalizations of the information matrix (IM) equality of White (1982).
The idea is that eigenspectrum-based statements of the IM equality reduce the degrees
of freedom of the test’s asymptotic distribution and lead to better size-power
properties, even in high dimensions. The gains are especially pronounced for
vine copulas, where additional benefits come from simplifications of score functions and the Hessian.
We derive the asymptotic distribution of the generalized
tests, accounting for the non-parametric estimation of the marginals and apply a 
parametric bootstrap procedure, valid when asymptotic critical values
are inaccurate. In Monte Carlo simulations, we study the behavior of the new
tests, compare them with several Cramer-von Mises type tests and confirm the
desired properties of the new tests in high dimensions.

The used R code for the Monte Carlo simulations in the Vine Copula part are given in this repo.
The used MATLAB code for the copula part can be requested from Yajing Zhu.

The code heavily builds on functionality from the **VineCopula package** ([CRAN](https://cran.r-project.org/web/packages/VineCopula/index.html), [github](https://github.com/ulf85/VineCopula)), in particular on the functions `RVineHessian` and `RVineGradient`.

## RVineGofTest_new.R

Since the code is not finally cleaned and brushed for a package or to be part of the VineCopula package the file naming and function names are quite lazily chosen. Sorry for that.  
Nevertheless, this file includes the main function `RVineGofTest_new2()`, which calculated the test statistics for various goodness-of-fit test developed in the paper. The function is inspired by the RVineGofTest function of the VineCopula package having also its argument names. The different goodness-of-fit (GOF) tests are referred by the `method` argument:

+ "White2" = Determinant White Test
+ "White3" = Trace White Test
+ "IR2" = Log Determinant IR Test
+ "log_trace" = Log Trace IMT
+ "log_GAIC" = Log Generalized Akaike Information Criterion (GAIC) IMT
+ "eigen" = Log Eigenspectrum IMT 
+ "eigen2" = Eigenvalue Test
 
Note: only a preliminary stage of the test statistics can be computed. 
Meaning that the necessary variance (-covariance matrix) has to be approximated
by bootstrapping and the test statistics have to be corrected.  
see correctTestStatistic()

## correctTestStatistic.R

The test statistics we get from `RVineGofTest_new2()` are only a preliminary 
stage of the test statistics. They have to be corrected with the variance
(covariance matrix).

In this file we have the function `correctTestStatistic()` with some sub-functions to do the "correction".
As input we assume a vector (or in some cases a matrix) of test statistics from a Bootstrap or Monte Carlo run.

## sizePower.R

Given the corrected test statistics we can approximate the p-value by Bootstrapping.
Given these (or asymptotic p-values) we can calculate the empirical size and power of the tests.  
We assume as input for the function `sizePower()` two vectors (or in some cases a matrices) of test statistics.
One for the true model (size) and one for the alternative (power).

## alternativeVineCopulaModels.R

In this file we give the code for the alternative Vine Copula models used in the paper for the Monte Carlo simulation. In the paper an R-vine is defined. The structure, the copula families and copula parameters are stored in an `RVM` object of the VineCopula package. Based on that `RVM` we can select and estimate the alternative Vine Copula models, namely a C-vine, a D-vine and an R-vine with just Gaussian copulas, i.e. a multivariate Gauss copula.
