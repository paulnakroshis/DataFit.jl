"""
    Δp(model, tᵢ, p, σₚ, α=1e-8)

Compute the partial derivatives of the model output at a specific time tᵢ, 
by numerically computing the partial derivative of the model with respect 
to the parameters p at time tᵢ using the central difference method.

Arguments:
- model: a function that takes an array of times and parameters 
         as input and returns an array of model outputs.
- tᵢ: the particular time at which to evaluate the model.
- p: a vector of parameters at which to evaluate the model.
- σₚ: a vector of parameter uncertainties from the bootstraped fits.
- α: a scalar specifying the fractional change in each parameter 
     used to compute the partial derivative, default is 1.0e-8. 
     Example: if  σ₁ = 0.1, then the step size, Δ, used to compute the central difference
     derivative with respect to σ₁ will be 0.1 * 1e-8.

Returns:
- ∂pᵢ: a vector containing the partial derivatives of the model 
       with respect to the model parameters at time tᵢ.

Usage:
- Calculate the numerical uncertainty of the model output at time t=2, 
  given a model function `my_model` and parameter values `p=[1,2,3]` 
  with uncertainties `σₚ=[0.1, 0.2, 0.3]`, using a step size of α=1.0e-6:
  ∂pᵢ = Δp(my_model, 2, p, σₚ, 1.0e-8)

Comments on method:
- the array  p₊ consists of an N x N array each row of which contains the fit parameters and 
  the diagonal element values each have a small added value δ=α*σₚ[i] where i is the row number. 
  Similarly for the array  p₋; except here, the diagonal elements have δ=α*σₚ[i] subtracted from the
  parameter value. 

"""
function Δp(model, tᵢ, p, σₚ, α)
    ∂pᵢ = zeros(length(p))clear
    ϵij = α*Diagonal(ones(length(p))) 
    p₊ = Matrix{Float64}(undef, length(p), length(p)) 
    p₋ = Matrix{Float64}(undef, length(p), length(p))
    for i in 1:length(p)
        p₊[i, : ] = p .+ ϵij[i, :].*σₚ[i] 
        p₋[i, : ] = p .-  ϵij[i, :].*σₚ[i]
        #use central difference method to estimate derivative
        ∂pᵢ[i] = (model(tᵢ, p₊[i, :]) .- model(tᵢ, p₋[i, : ]))/(2*α.*σₚ[i])
    end
    return ∂pᵢ
end



"""
    uncertainty(model, t::Vector, σₜ::Vector, p::Vector, σₚ::Vector, α=1e-6)

Computes the uncertainty in a model fit to a time-series data set.
Arguments:
- model: a function that takes an array of times and parameters 
         as input and returns an array of model outputs.
- t :     a vector of time values at which to evaluate the model.
- σₜ :     a vector of time uncertainties 
- p :     a vector of parameters at which to evaluate the model.
- σₚ :    a vector of parameter uncertainties from the bootstraped fits.
- α :     a scalar specifying the fractional change in each parameter 
         used to compute the partial derivative, default is 1.0e-8. 
         Example: if  σ₁ = 0.1, then the step size, Δ, used to compute the central difference
         derivative with respect to σ₁ will be 0.1 * 1e-8.

Returns:
- Δy:    a vector containing the uncertainty of the model fit at each time.
         The error is computed by adding the errors in quadrature. The unique feature of 
         this is that this uncertainty includes the effects of both time and parameter
         uncertainty. 
"""
function uncertainty(model, t::Vector, σₜ::Vector, p::Vector, σₚ::Vector, α=1e-6)
    Δy = zeros(length(t))
    for i in 1:length(t)
        if i==1
            ∂ₜ = (model(t[i] + α*σₜ[i], p) - model(t[i],p))/(α*σₜ[i])
            ∂ₚ = Δp(model, t[i], p, σₚ, α)
            Δy[i] = sqrt((∂ₜ*σₜ[i])^2 + sum((∂ₚ .* σₚ).^2))
        elseif i<length(t)
            ∂ₜ = (model(t[i] + α*σₜ[i], p) - model(t[i]-α*σₜ[i],p))/(2*α*σₜ[i])
            ∂ₚ = Δp(model, t[i], p, σₚ, α)
            Δy[i] = sqrt((∂ₜ*σₜ[i])^2 + sum((∂ₚ .* σₚ).^2))
        else
            ∂ₜ = (model(t[i], p) - model(t[i]-α*σₜ[i],p))/(α*σₜ[i])
            ∂ₚ = Δp(model, t[i], p, σₚ, α)
            Δy[i] = sqrt((∂ₜ*σₜ[i])^2 + sum((∂ₚ .* σₚ).^2))
        end
    end
    return Δy
end