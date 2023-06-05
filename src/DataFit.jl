module DataFit

include("dataplot.jl")
export "plot_defaults", "plot_raw_data"

include("fit.jl")
export fit_data,  spacing_info, bootstrap_fit
export  Δp, uncertainty

end