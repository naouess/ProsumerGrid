using JLD2, FileIO, GraphIO, CSV, DataFrames

begin
	## Generating demand function out of real minigrid data
	dem_data_week = CSV.read("$dir/demand/Week_March02_March08_2020.csv")

	load = t->dem_data_week[!,:load][Int(floor(mod(t,24*3600*7) / 900)+1)]
	pv_infeed = t->dem_data_week[!,:pv_infeed][Int(floor(mod(t,24*3600*7) / 900)+1)]

	periodic_demand = t -> vcat(load(t), pv_infeed(t)) ./ 30. # , week_G4(t), (week_winter(t) + week_G1(t) + week_G4(t))./3)

	# For additional fluctuations:
	# samples = 24*60
	# inter = interpolate([10. ./100 * randn(2) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	# residual_demand = t -> inter(1. + t / (24*3600) * samples)

	demand = t -> periodic_demand(t) # + residual_demand(t)

end
