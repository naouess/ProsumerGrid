using JLD2, FileIO, GraphIO, CSV, DataFrames

begin
	struct demand_amp_var
		demand
	end

	function (dav::demand_amp_var)(t)
		index = Int(floor(t / (24*3600)))
		dav.demand[index + 1,:]
	end

	## Generating demand function out of Standard load profile H0
	dem_data = CSV.read("$dir/demand/profil.csv")
	dem_data_week = CSV.read("$dir/demand/profil_week.csv")

	week_winter = t->dem_data_week[!,:Winterwoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
	week_summer = t->dem_data_week[!,:Sommerwoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
	week_between = t->dem_data_week[!,:Uebergangswoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
	week_G1 = t->dem_data_week[!,:WinterwocheG1][Int(floor(mod(t,24*3600*7) / 900)+1)]
	week_G4 = t->dem_data_week[!,:WinterwocheG4][Int(floor(mod(t,24*3600*7) / 900)+1)]

	week = t->vcat(week_winter(t), week_G1(t), week_G4(t), (week_winter(t) + week_G1(t) + week_G4(t))./3)


	periodic_demand = t ->  week(t) ./100
	samples = 24*60
	inter = interpolate([10. ./100 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	residual_demand = t -> inter(1. + t / (24*3600) * samples)

end
