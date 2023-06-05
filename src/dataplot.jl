"""
   plot_defaults(w)
   
Sets nice LaTeX-like plot default values. 
For reference:  
a4 paper is 210 x 297 mm  
or 2480 x 3508 pixels at 300 dpi 
"""
function plot_defaults(width, height) # plot width and height in pixels
	default(fontfamily = "Computer Modern", 
			size=(width, height),
			titlefont = (16), 
			legendfontsize = 10, 
			guidefont = (16, :black), 
			tickfont = (12, :black), 
	        framestyle = :box, 
			yminorgrid = true, 
			xminorgrid= true,
			legend = :outertopright, 
			dpi=300, 
			margin = 7mm)
end

"""
plot_raw_data(tdata, ydata, 
   				xerr, yerr, 
				xlabel, ylabel, 
				legend_label, 
				smoothing=1)

	Plot raw data with error bars and optional smoothing.

	Arguments:
	- `tdata`: time data
	- `ydata`: y data
	- `xerr`: x error bars
	- `yerr`: y error bars
	- `xlabel`: x axis label
	- `ylabel`: y axis label
	- `legend_label`: legend label
	- `smoothing`: number of passes of binomial smoothing to apply to data. Default is 1.

	Returns:
	- `p`: plot object
"""
function plot_raw_data(tdata, ydata, xerr, yerr, 
						xlabel, ylabel, legend_label, smoothing=1)
	
	if typeof(smoothing)≠Int64 || smoothing < 0
		error("smoothing must be a positive integer")
	end
	if smooting	> 0 				
		p = scatter(tdata, ydata, 
					xerr=xerr, yerr=yerr, 
					xlabel = L"xlabel", ylabel =L"ylabel",
					marker=:circle, markersize=2, markeralpha=0.35, markercolor=:brown,
					label = L"legend_label",  
					legend = :outertopright, 
					dpi=300,
					margin = 7mm
					)
			plot!(tdata,  Smoothing.binomial(ydata, smoothing), 
				 alpha=0.6, color=:darkred, label="$smoothing pass binomial smoothing")
		
	elseif smoothing == 0
		p = scatter(tdata, ydata, 
					xerr=xerr, yerr=yerr, 
					xlabel = L"xlabel", ylabel =L"ylabel",
					marker=:circle, markersize=2, markeralpha=0.35, markercolor=:brown,
					label = L"legend_label",  
					legend = :outertopright, 
					dpi=300,
					margin = 7mm
					)
		return p
	end
end
