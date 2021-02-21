using JLD2, FileIO, GraphIO, CSV, DataFrames
using Statistics
using StatsBase

## 02.03 + 03.03.2020 and 5 Diesel days ?? which data did I use for day 1-5 ?
begin
	## Generating demand function out of real minigrid data
	dem_data_twoweeks = CSV.read("$dir/src/demand/TwoWeeks.csv")

	load = t->dem_data_twoweeks[!,:load][Int(floor(mod(t,24*3600*num_days) / 900)+1)]
	pv_infeed = t->dem_data_twoweeks[!,:pv_infeed][Int(floor(mod(t,24*3600*num_days) / 900)+1)]
	diesel = t -> dem_data_twoweeks[!,:diesel][Int(floor(mod(t,24*3600*num_days) / 900)+1)]
	demand_real = t -> vcat(load(t), pv_infeed(t)) ./ 50.
end

plot(0:num_days*l_day, t -> diesel(t) ./ 50)

for (i,t) in enumerate(sol.t) #sol.t #1:24*num_days+1

	hour = Int(div(t, 3600) + 1)

	ILC_power_node[i] = ILC_node[hour]
	LI_power_node[i] = - KP[node] * sol(t)[indices_ω[node]] + sol(t)[indices_χ[node]]

	control_node[i] = LI_power_node[i] + ILC_power_node[i]

	global control_integral_node += abs(control_node[i]) * delta_t
	# global control_solt_integral += abs(control_solt[i])
end

using Plots
begin
	plot(120*3600:3600:num_days*l_day, t -> diesel(t) ./ 50, labels="Data", linecolor=:lightblue,
	fill=(0, 0, :lightblue),
	# alpha=1,
	#linestyle=:dot,
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm,
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	# ribbon=mismatch_solt
	)
	ylims!(-1,1)
	plot!(120*3600:3600:24*num_days*3600, LI_exact[120:num_days*24,3] .+ vcat(ILC_power[:,:,3]'...)[120:num_days*24],
	#mismatch_solt,
	fill=(0, 0, :black),
	# xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	linecolor=:black, linestyle=:solid,
	#linewidth=1
	)
	savefig("$dir/plots/VII_diesel.png")
end

begin
	plot(1:3600:24*num_days*3600, LI_exact[1:num_days*24,3] .+ vcat(ILC_power[:,:,3]'...)[1:num_days*24],
	#mismatch_solt,
	fill=(0, 0, :black),
	# xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	linecolor=:black, linestyle=:solid,
	# alpha=1,
	#linestyle=:dot,
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm,
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	# ribbon=mismatch_solt
	)
	ylims!(-1,1)
	# savefig("$dir/plots/VII_diesel_simulated.png")
end
# begin
# 	struct demand_amp_var
# 		demand
# 	end
#
# 	function (dav::demand_amp_var)(t)
# 		index = Int(floor(t / (24*3600)))
# 		dav.demand[index + 1,:]
# 	end
#
#     ## Generating synthetic demand function:
# 	# slowly increasing and decreasing amplitude - only working for <= 20 days now
# 	demand_amp1 = demand_amp_var(repeat([80 80 80 10 10 10 40 40 40 40 40], outer=Int(N/4))') # random positive amp over days by 10%
# 	demand_amp2 = demand_amp_var(repeat([10 10 10 80 80 80 40 40 40 40 40], outer=Int(N/4))') # random positive amp over days by 10%
# 	demand_amp3 = demand_amp_var(repeat([60 60 60 60 10 10 10 40 40 40 40], outer=Int(N/4))') # random positive amp over days by 10%
# 	demand_amp4 = demand_amp_var(repeat([30 30 30 30 10 10 10 80 80 80 80], outer=Int(N/4))') # random positive amp over days by 10%
# 	demand_amp = t -> vcat(demand_amp1(t), demand_amp2(t), demand_amp3(t), demand_amp4(t))
#
# 	periodic_demand =  t-> demand_amp(t)./100 .* sin(t*pi/(24*3600))^2
# 	samples = 24*4
# 	inter = interpolate([.2 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
# 	residual_demand = t -> inter(1. + t / (24*3600) * samples)
# 	demand_synth = t -> periodic_demand(t) + residual_demand(t)
# end
