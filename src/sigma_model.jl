# numerical uncertainty for the model at time tᵢ computed using the partial derivative for the model with respect to the parameters. p,  at time tᵢ, and adding errors in quadrature.
function δy_at_tᵢ(model, tᵢ, p, σₚ, α=0.001)
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