
<!-- README.md is generated from README.Rmd. Please edit that file -->

# arima2

<!-- badges: start -->
<!-- badges: end -->

The goal of `arima2` is to provide a set of tools to aid in the analysis
of time series data in `R`. One such function within the `arima2`
package is called `arima2`, which provides an interface to estimating
Auto Regressive Integrated Moving Average (ARIMA) models using a
random-restart algorithm. This function improves on the functionality of
the `stats::arima` method, as it has the potential to increase the
likelihood of the final output model. By design, the function cannot
result in models with lower likelihoods than that of the `stats::arima`
function. The potential for increased model likelihoods is obtained at
the cost of computational efficiency, as the function is order $O(n)$
times slower than the `stats::arima` function, where $n$ is the number
of random restarts. Because the estimation of ARIMA models takes only a
fraction of a second, we are of the opinion that potential to increase
model likelihoods is well worth this computational cost. The `arima2`
function is implemented by modifying the source code of the
`stats::arima` function.

## TODO:

- [x] Implement `arima2`
- [ ] Create simulate function for ARIMA models.
- [ ] Add datasets to the package
- [ ] (Dhajanae) `ggplot2` figures for ARIMA objects
- [ ] (Noel) Create function that makes AIC table.
- [ ] `probe`: a function to compare simulated data to a model
- [ ] (Noel)`auto.arima2:` implements that `auto.arima` function using
  `arima2`.
- [ ] (Jesse) `polyroots:` A function to get the roots of the
  polynomials created by the `AR` and `MA` coefficients of a model.

## Installation

You can install the development version of arima2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jeswheel/arima2")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(arima2)
## basic example code
```

Don’t forget to render the `README.Rmd` if you update it. If you don’t,
then there is a GitHub Action implemented to re-render the `README.Rmd`,
but this is less ideal.
