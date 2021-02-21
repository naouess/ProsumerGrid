using JLD2, FileIO, GraphIO, CSV, DataFrames
# using Statistics
# using StatsBase

begin
	## Generating demand function out of real minigrid data
	dem_data_twoweeks = CSV.read("$dir/src/data/TwoWeeks.csv", DataFrame)

	load = t->dem_data_twoweeks[!,:load][Int(floor(mod(t,24*3600*num_days) / 900)+1)]
	pv_infeed = t->dem_data_twoweeks[!,:pv_infeed][Int(floor(mod(t,24*3600*num_days) / 900)+1)]

	demand_real = t -> vcat(load(t), pv_infeed(t)) ./ 50.
end
