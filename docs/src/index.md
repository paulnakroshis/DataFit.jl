```@meta
CurrentModule = DataFit
```

# DataFit

Documentation for [DataFit](https://github.com/paulnakroshis/DataFit.jl).

```@index
## Background/Introduction
This package is designed to enable easy fitting of time series data  to a user-supplied model 
and to take into account the experimental uncertainties (in both the independent and dependent 
variables) in the estimation of model parameter uncertainties. 

One area this package could be helpful is in undergraduate laboratories where students are 
reminded of the importance of including uncertainties, but when fitting their data to a model,
they do not know how to include the uncertainties on their measured data in performing their
model fit. 

## Package features
The user must supply a five items:  

    1. a set of time measurements as a vector
    2. the uncertainties on these time measurements as a vector 
    3. a set of dependent variable measurements
    4. the uncertainties on the dependent variable measurements
    5. a model function used to fit the data
    6. the number of bootstrapped data sets to then create and fit

DataFit.jl uses the CurveFit.jl package to obtain the best fit to the data.
The extra feature of the this package is that the user then specifies how many bootstrapped
data sets to create to then provide separate estimates for the best fit parameters and their 
uncertainty. 

The package also provides a simple publication ready plot using $\LaTeX$ axes labels and 
computer modern fonts; this plot includes:  

    a. the raw data with error bars
    b. the best fit
    c. a 1 sigma (or 2,3,4,or 5 sigma) error band 
    d. a residuals plot above the main plot

### Bootstrapped data
A bootstrapped data set is a data set of the same size as the original, with points randomly chosen from 
a Gaussian error ellipsoid centered on the the orginal data points and whose standard deviations, $\sigma_t$ 
and $\sigma_y$ (using $y$ to represent the dependent variable) are given by the experimental uncertainties
in the independent and dependent variables. There is an option to use a multiple of $\sigma_t$ 
and $\sigma_y$ in this bootstrapping process. 

Once this bootstrapped data set is created, a fit is done to the user's model function, and this process is repeated
M times; the routine then returns the mean and the standard deviation of the best fit parameters as a new estimate of 
model parameters which thus takes into account the experimental uncertainties. 

Of course in the limit of infinite bootstrapped samples, the average parameter values will equal the best fit parameter
values, so really the benefit of the bootstrap method is to give more realistic estimates of the uncertainties on the 
fit parameters. Nonetheless, the bootstrap routine returns both the mean and the standard deviation of the model parameters.  

### Computing fit uncertainties
In order to **plot** the effect of the model parameter uncertainties as determined by the bootstrapped method, 
we need to know how the uncertainties effect the model fit value. To do so we need to only evaluate the partial derivative of 
the model with respect to time and each model parameter at each data point being considered. Once we have these partial 
derivatives, the model uncertainty, $\Delta y$, is computed by adding the errors in quadrature:  

$\Delta y(t_i) =  \sqrt{\left(\frac{\partial y}{\partial t}\sigma_t\right)^2 + 
   \sum_{j=1}^N \left(\frac{\partial y}{\partial p_j} \sigma_{p_j}\right)^2},$

where the number of model parameters is $N$, and $\sigma_t$ and $\sigma_{p_j}$ are the arrays containing 
uncertainties in the times and parameters values. 

The array $\Delta y$ can then be used to make a **ribbon plot** (an option in Plots.jl). 


```

```@autodocs
Modules = [DataFit]
```
