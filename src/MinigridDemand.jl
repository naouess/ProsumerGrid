using JLD2, FileIO, GraphIO, CSV, DataFrames

begin
	## Generating demand function out of real minigrid data
	dem_data_week = CSV.read("$dir/demand/profil_week.csv")

	load = t->dem_data_week[!,:load][Int(floor(mod(t,24*3600*7) / 900)+1)]
	pv_infeed = t->dem_data_week[!,:pv_infeed][Int(floor(mod(t,24*3600*7) / 900)+1)]

	demand = t -> vcat(load(t), pv_infeed(t)) ./ 5. # , week_G4(t), (week_winter(t) + week_G1(t) + week_G4(t))./3)

	# For additional fluctuations:
	# samples = 24*60
	# inter = interpolate([10. ./100 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	# residual_demand = t -> inter(1. + t / (24*3600) * samples)

end
