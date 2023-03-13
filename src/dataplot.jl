"""
   plot_defaults(w)
   
Sets nice LaTeX-like plot default values. 
"""
function plot_defaults()
	default(fontfamily = "Computer Modern", size=(1200,800), titlefont = (16), legendfontsize = 10, 
	        guidefont = (16, :darkgreen), tickfont = (12, :black), 
	        framestyle = :box, yminorgrid = true, xminorgrid= true,legend = :outertopright, dpi=600, margin = 7mm)
