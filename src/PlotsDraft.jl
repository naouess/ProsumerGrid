begin
# Plotting the demand
	load_amp_hourly = [maximum.(dd(t)) for t in 1:3600:3600*24*num_days]
	plot(0:7*l_day, t -> dd(t)[1],ytickfontsize=14,
	               xtickfontsize=18, margin=5Plots.mm,
	    		   legendfontsize=12, linewidth=3,xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)),xaxis=("days [c]",font(14)), yaxis=("normed demand",font(14)), legend=nothing)
end
begin
# hourly plotting
	plot(1:num_days*24, ILC_power_hourly[1:24*num_days] , legend=nothing, label=L"$ u_j^{ILC}$", ytickfontsize=14,
	               xtickfontsize=18,
	    		   legendfontsize=12, linewidth=3,xaxis=("time [h]",font(14)), yaxis=("normed power",font(14)))
	plot!(1:num_days*24+1, mean(hourly_energy, dims=2)/3600 , label=L"y^{c,h}", linewidth=3)
	# plot!(1:24*num_days, mean.(load_amp_hourly), label = "peak demand", linewidth=3)
end

# second-wise
begin
	plot(1:3600:num_days*24*3600,  ILC_power_hourly[1:num_days*24]./ maximum(ILC_power_hourly), label=L"$P_{ILC, j}$", ytickfontsize=14,
	               xtickfontsize=18,
	    		   legendfontsize=10, linewidth=3,xaxis=("time [s]",font(14)), yaxis=("normed quantities [a.u.]",font(14)))
	plot!(1:3600:24*num_days*3600,mean(hourly_energy[1:num_days*24], dims=2) ./ maximum(hourly_energy), label=L"y_h",linewidth=3, linestyle=:dash)
	# plot!(0:num_days*l_day, t -> dd(t)[1], label = "demand",linewidth=3, alpha=0.3)
end

# To reproduce figure 6 in paper:
begin
	load_amp_hourly_N = [dd(t) for t in 1:3600:3600*24*num_days]
	load_amp_hourly = sum.(load_amp_hourly_N)
	load_amp_daily = sum.(mean(reshape(load_amp_hourly, 24,num_days)',dims=2))
end
begin
	plot(2:num_days, sum(ILC_power_agg[2:num_days,1,:],dims=2), label=L"$\bar u^{ILC}$", ytickfontsize=14,
	               xtickfontsize=18, margin=5Plots.mm,
	    		   legendfontsize=14, linewidth=3,xaxis=("days [c]",font(14)), yaxis=("normed power",font(14)), legend=:right)
	plot!(2:num_days, sum(mean_energy_d[2:num_days],dims=2) ./ 3600, label=L"\bar y^{c}", linewidth=3, linestyle=:dash)
	plot!(2:num_days, load_amp_daily[2:num_days] , label = L"\bar P^d", linewidth=3, linestyle=:dashdot)
end
begin
	plot(2:num_days, sum(ILC_power_agg_norm[2:num_days,1,:],dims=2), label=L"$\bar u^{ILC}$", ytickfontsize=14,
	               xtickfontsize=18, margin=5Plots.mm,
	    		   legendfontsize=14, linewidth=3,xaxis=("days [c]",font(14)), yaxis=("normed power",font(14)), legend=:right)
	plot!(2:num_days, sum(norm_energy_d[2:num_days],dims=2) ./ 3600, label=L"\bar y^{c}", linewidth=3, linestyle=:dash)
	plot!(2:num_days, load_amp_daily[2:num_days] , label = L"\bar P^d", linewidth=3, linestyle=:dashdot)
end
