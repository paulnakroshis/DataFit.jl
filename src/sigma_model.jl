"""
δy_at_tᵢ(model, tᵢ, p, σₚ, α=1.0e-6)

Compute the numerical uncertainty of the model output at time tᵢ, 
by computing the partial derivative of the model with respect 
to the parameters p at time tᵢ and adding the errors in quadrature.

Arguments:
- model: a function that takes an array of times and parameters as input and returns an array of model outputs.
- tᵢ: the particular time at which to evaluate the model.
- p: a vector of parameters at which to evaluate the model.
- σₚ: a vector of parameter uncertainties from the bootstraped fits.
- α: a scalar specifying the fractional change in each parameter 
     used to compute the partial derivative, default is 1.0e-6. 
     Example: if  σ₁ = 0.1, then the step size, Δ, used to compute the central difference
     derivative with respect to σ₁ will be 0.1 * 1e-6.

Returns:
- dy: a scalar representing the numerical uncertainty of the model output at time tᵢ.

Examples:
- Calculate the numerical uncertainty of the model output at time t=2, 
given a model function `my_model` and parameter values `p=[1,2,3]` 
    with uncertainties `σₚ=[0.1, 0.2, 0.3]`, using a step size of α=1.0e-6:
  δy = δy_at_tᵢ(my_model, 2, p, σₚ, 1.0e-6)




"""
function δy_at_tᵢ(model, tᵢ, p, σₚ, α=1.0e-6)
    Δ = α.*σₚ
    ∂pᵢ = zeros(length(p))
     p₊ = p .+ Δ
     p₋ = p .-  Δ
    ∂pᵢ = (model([tᵢ], p₊) .- model([tᵢ], p₋))/(2*Δ)
    dy = sqrt(sum(∂pᵢ.^2 .* σₚ.^2))
    return dy
end

# uncertainty in the model function at each time point
function dy(model, t, p, σₚ; α=0.001)
    i=1
    Δy = zeros(length(t))
    for  tᵢ in t    
        Δy[i] = δy_at_tᵢ(model, tᵢ, p, σₚ, α)
        i += 1
    end
    return  Δy
end

