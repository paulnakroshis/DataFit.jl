# DataFit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://paulnakroshis.github.io/DataFit.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://paulnakroshis.github.io/DataFit.jl/dev/)
[![Build Status](https://github.com/paulnakroshis/DataFit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/paulnakroshis/DataFit.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/paulnakroshis/DataFit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/paulnakroshis/DataFit.jl)

The goal of this is to create a Julia package which will allow the user to
fit a data set with one independent variable to a model with N parameters, and
to do so while taking accound the uncertainties in the independent and dependent
variables.

I'll use a time series as an example to frame the discussion.

The user supplies

 1. a time sequence together with the uncertainty at each time,
 2. a set of measurements, [y] together with their uncertainties,
 3. a model function with N parameters with which to fit the data
 4. whether to smooth the data before fitting
the code will return the best fit model parameters, and will determine the
uncertainty by creating M bootstrapped data sets and refitting to each. The
uncertainty in the N parameters is determined by the standard deviation of
the M bootstrapped data set parameter values.

Then, this code is used to plot the data along with 1$\sigma$ and 3$\sigma$
confidence bands (not sure if this is the right term for this).
