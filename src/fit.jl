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
An error will be thrown if the vector is not strictly increasing. The reason for this test is that the
procedure used in `bootstrap_fit` requires that the data be well-spaced; if the user creates a bootstrapped
data set that violates this requirement, I want to know about it.

# Arguments
- `v::Vector`: Vector of data to be analyzed.
- `full_stats::Bool`: If `true`, returns the minimum, maximum, mean, 
                    and standard deviation of the spacing between elements of `v`. 
                    The default is `false`, and then the function only returns the minimum spacing.

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
ERROR: Vector must be strictly increasing
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
        error("Vector must be strictly increasing")
    end
    Δtmin = minimum(Δt)
    Δtmax = maximum(Δt)
    Δtmean = mean(Δt)
    Δtstd = std(Δt)
    if full_stats
        Δtmin = minimum(Δt)
        Δtmax = maximum(Δt)
        Δtmean = mean(Δt)
        Δtstd = std(Δt)
        return round.([Δtmin, Δtmax, Δtmean, Δtstd]; digits=6)
    else
        return round(Δtmin; digits=6)
    end
end



function bootstrap_fit(model, tdata, δt, Σt, ydata, δy, Σy, p0, m_samples)
    fit_params, σ = fit_data(model, tdata, ydata, p0)
    fit_params_array = zeros(length(fit_params), m_samples)
    
    for i in 1:m_samples
        t_boot = rand.(Normal.(tdata, abs.(Σt*δt))) # create bootstrap time data
        y_boot = rand.(Normal.(ydata, abs.(Σy*δy))) # create bootstrap y data
        fit_params_array[:,i] = fit_data(model, t_boot, y_boot, p0)[1]
    end
    p = vec(mean(fit_params_array, dims=2))
    σₚ =vec(std(fit_params_array, dims=2))
    return p, σₚ
end