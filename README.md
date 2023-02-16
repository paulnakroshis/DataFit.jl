# DataFit.jl
## 16 Feb 2023: Not ready for use yet
### In the docs folder is a notebook used for development.

The goal of this Julia package is to allow the user to fit a data set with one independent variable to a model with N parameters, and to do so while taking accound the uncertainties in the independent and dependent 
variables. 

I'll use a time series as an example to frame the discussion. 

The user supplies 

     1. a time sequence together with the $1\sigma$ uncertainties at each time,
     2. a set of measurements, [y] together with their uncertainties,
     3. a model function with N parameters with which to fit the data,
     4. the number of bootstrapped sample data sets to create,
     5. whether to use $1\sigma$, $2\sigma$, ... or more errors in y and t when bootstrapping.

the code will return the best fit model parameters, and will determine the 
uncertainty by creating M bootstrapped data sets and refitting to each. The 
uncertainty in the N parameters is determined by the standard deviation of 
the M bootstrapped data set parameter values. 

Then, this code is used to plot the data along with error bars and confidence band(s)
