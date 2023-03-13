"""
fit_data(model, tdata, ydata, p0)

Fit the model function `model` to the data `ydata` at time points `tdata`
using initial parameters `p0`.

# Arguments
- `model::Function`: The model function to fit.
- `tdata::Vector`: A vector of time points.
- `ydata::Vector`: A vector of corresponding data points.
- `p0::Vector`: A vector of initial parameter values.

# Returns
A tuple `(fit_params, σ)` of the fitted parameter values and their standard deviations.

# Example
```julia
model(x, p...) = p[1] .* exp.(-p[2] .* x) .* sin.(p[3] .* x .+ p[4])
tdata = range(0, 10, length=100)
ydata = model(tdata, 1.0, 0.5, 2π, π/2) + randn(size(tdata))
p0 = [1.0, 1.0, 2π, π/2]
fit_params, σ = fit_data(model, tdata, ydata, p0)
"""
function fit_data(model, tdata, ydata, p0)
    fit = curve_fit(model, tdata, ydata, p0)
    fit_params = fit.param
    σ = stderror(fit)
    return fit_params, σ # returns a tuple of fit parameters and their standard deviations
end


"""
    spacing_info(v::Vector, full_stats::Bool=false)

Returns the minimum spacing between elements of a vector `v`. If `full_stats` is `true`, 
returns the minimum, maximum, mean, and standard deviation of the spacing between elements of `v`.
The default is `false`, and then the function only returns the minimum spacing.

An error will be thrown if the vector is not strictly increasing. The reason for this test is that 
the procedure used in `bootstrap_fit` requires that the data be time ordered; if the user creates a 
bootstrapped data set that violates this requirement, I want to know about it.

# Arguments
- `v::Vector`: Vector of data to be analyzed.
- `full_stats::Bool`: If `true`, returns the minimum, maximum, mean, 
                    and standard deviation of the spacing between elements of `v`. 
                    The default is `false`, and in this case the function only returns the                         minimum spacing.
---
# Examples
## Basic usage with well-spaced data:
```julia-repl
julia> using Statistics
julia> include("fit.jl")
julia> t = [1.0, 2.0, 3.0, 3.5, 4.2, 5.0, 6.1, 7.4, 7.6];

julia> spacing_info(t)
0.2

julia> spacing_info(t, true)
4-element Vector{Float64}:
 0.2
 1.3
 0.825
 0.353553
```

## Basic usage with poorly-spaced data:
```julia-repl
julia> t = [1.0, 2.0, 3.0, 3.5, 3.4, 5.0, 6.1, 7.4, 7.6];

julia> spacing_info(t)
ERROR: Vector must be strictly increasing. Your generated 
bootstrap array is not time ordered. Try reducing the number 
of standard deviations used in your gaussian sample. For regularly 
spaced time samples the value of nΣt should be 1.0; 
did you override this?
Stacktrace:
 [1] error(s::String)
   @ Base ./error.jl:35
 [2] spacing_info(v::Vector{Float64}, full_stats::Bool)
   @ Main ~/Dropbox/DocumentsF/Github/DataFit.jl/src/fit.jl:32
 [3] spacing_info(v::Vector{Float64})
   @ Main ~/Dropbox/DocumentsF/Github/DataFit.jl/src/fit.jl:30
 [4] top-level scope
   @ REPL[41]:1
```
"""
# analyze spacing bewteen elements of a Vector
function spacing_info(v::Vector, full_stats::Bool=false)
    Δt =v[2:end] - v[1:end-1]
    if any(Δt .<= 0)
        error("Vector must be strictly increasing. 
              Your generated bootstrap array is not time ordered. 
              Try reducing the number of standard deviations used in 
              your gaussian sample. For regularly spaced time samples
              the value of nΣt should be 1.0; did you override this?")
    end
    
    Δtmin = minimum(Δt)
    
    if full_stats
        Δtmax = maximum(Δt)
        Δtmean = mean(Δt)
        Δtstd = std(Δt)
        return round.([Δtmin, Δtmax, Δtmean, Δtstd]; digits=6)
    else
        return round(Δtmin; digits=6)
    end
end

"""
bootstrap_fit(model, tdata, δt, Σt, ydata, δy, Σy, p0, m_samples)

Bootstrap fitting of data to a model, returning the mean and standard deviation of the fit parameters.

# Arguments
- `model`: A function that takes a time vector and a parameter vector as input and returns a predicted output vector.
- `tdata`: A vector of the observed time points.
- `δt`: A vector of the uncertainties in the observed time points.
- `ydata`: A vector of the observed output values.
- `δy`: A vector of the uncertainties in the observed output values.
- `p0`: A vector of initial parameter values for the model.
- `m_samples`: The number of bootstrap samples to generate.
- `nΣt`: A scalar representing what multiple of δt you want to use for the width of the Gaussian used to generate fake t data. This should be set to 1.0 for regularly spaced data-- otherwise the time series will not be time ordered. For irregularly spaced data, it is possible to choose values larger than 1.0. In any case, there is a check built into the function to warn the user if the time series data is unordered.
- `nΣy`: A scalar representing what multiple of δy you want to use for the width of the Gaussian used to generate fake y data. Default value=1.0.

# Returns
- `p`: A vector of the mean values of the fit parameters across all bootstrap samples.
- `σₚ`: A vector of the standard deviations of the fit parameters across all bootstrap samples.
"""
function bootstrap_fit(model, tdata, δt, ydata, δy, p0, m_samples, nΣt=1.0, nΣy=1.0)
    fit_params, σ = fit_data(model, tdata, ydata, p0)
    fit_params_array = zeros(length(fit_params), m_samples)
    
    for i in 1:m_samples
        t_boot = rand.(Normal.(tdata, abs.(Σt*δt))) # create bootstrap time data
        y_boot = rand.(Normal.(ydata, abs.(Σy*δy))) # create bootstrap y data
        fit_params_array[:,i] = fit_data(model, t_boot, y_boot, p0)[1]
    end
    p   = vec(mean(fit_params_array, dims=2))
    σₚ = vec(std(fit_params_array, dims=2))
    return p, σₚ
end