function fit_data(model, tdata, ydata, p0)
    fit = curve_fit(model, tdata, ydata, p0)
    fit_params = fit.param
    σ = stderror(fit)
    return fit_params, σ # returns a tuple of fit parameters and their standard deviations
end


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
        return Δtmin, Δtmax, Δtmean, Δtstd
    else
        return Δtmin
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