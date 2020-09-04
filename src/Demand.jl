using JLD2, FileIO, GraphIO, CSV, DataFrames

begin
	struct demand_amp_var
		demand
	end

	function (dav::demand_amp_var)(t)
		index = Int(floor(t / (24*3600)))
		dav.demand[index + 1,:]
	end

	## Generating demand function as in sim_demand_learning.jl

	# dem_data = CSV.read("$dir/profil.csv")
	#
	# dem_data_week = CSV.read("$dir/profil_week.csv")
	#
	# week_winter = t->dem_data_week[!,:Winterwoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
	# week_summer = t->dem_data_week[!,:Sommerwoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
	# week_between = t->dem_data_week[!,:Uebergangswoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
	# week_G1 = t->dem_data_week[!,:WinterwocheG1][Int(floor(mod(t,24*3600*7) / 900)+1)]
	# week_G4 = t->dem_data_week[!,:WinterwocheG4][Int(floor(mod(t,24*3600*7) / 900)+1)]
	#
	# week = t->vcat(week_winter(t), week_G1(t), week_G4(t), (week_winter(t) + week_G1(t) + week_G4(t))./3)
	#
	#
	# periodic_demand = t ->  week(t) ./100
	# samples = 24*60 #24*4
	# inter = interpolate([10. ./100 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	# residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

    ## Generating synthetic demand function:

	# slowly increasing and decreasing amplitude - only working for <= 20 days now
	demand_amp1 = demand_amp_var(repeat([80 80 80 10 10 10 40 40 40 40 40], outer=Int(N/4))') # random positive amp over days by 10%
	demand_amp2 = demand_amp_var(repeat([10 10 10 80 80 80 40 40 40 40 40], outer=Int(N/4))') # random positive amp over days by 10%
	demand_amp3 = demand_amp_var(repeat([60 60 60 60 10 10 10 40 40 40 40], outer=Int(N/4))') # random positive amp over days by 10%
	demand_amp4 = demand_amp_var(repeat([30 30 30 30 10 10 10 80 80 80 80], outer=Int(N/4))') # random positive amp over days by 10%
	demand_amp = t -> vcat(demand_amp1(t), demand_amp2(t), demand_amp3(t), demand_amp4(t))

	periodic_demand =  t-> demand_amp(t)./100 .* sin(t*pi/(24*3600))^2
	samples = 24*N
	inter = interpolate([.2 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range
end
# plot(1:24*3600, t -> periodic_demand(t)[1] + residual_demand(t)[1])
