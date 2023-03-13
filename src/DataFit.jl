module DataFit

export fit_data,  spacing_info, bootstrap_fit
include("fit.jl")

export  Δp, uncertainty
include("sigma_model.jl")



end
