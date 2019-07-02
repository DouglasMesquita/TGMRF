
# TGMRF

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/FLAMES)](https://cran.r-project.org/package=FLAMES) -->

<!-- [![Travis build status](https://travis-ci.org/DouglasMesquita/FLAMES.svg?branch=master)](https://travis-ci.org/DouglasMesquita/FLAMES) -->

<!-- [![Codecov test coverage](https://codecov.io/gh/DouglasMesquita/FLAMES/branch/master/graph/badge.svg)](https://codecov.io/gh/DouglasMesquita/FLAMES?branch=master) -->

<!-- badges: end -->

## Overview

TGMRF: Transformed Gaussian Markov Random Fields is a package that allow
the use of several approaches to fit Poisson intensities.

## Installation

``` r
# Install from CRAN (when available)
install.packages("TGMRF")
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("DouglasMesquita/TGMRF")
```

## Usage

`library(TGMRF)` will load **tgmrf**, a function to fit Poisson
regression models. The user can choose for several models to model the
poisson intensities: **gamma-var**, **gamma-shape**, **gamma-scale** and
**log-normal**.

See `?tgmrf` for examples.