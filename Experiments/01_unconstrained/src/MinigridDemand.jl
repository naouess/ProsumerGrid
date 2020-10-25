using JLD2, FileIO, GraphIO, CSV, DataFrames
using Statistics
using Random
Random.seed!(42)
using Interpolations

## Real demand
begin
	## Generating demand function out of real minigrid data
	dem_data_twoweeks = CSV.read("$dir/src/demand/TwoWeeks.csv")

	load = t -> dem_data_twoweeks[!,:load][Int(floor(mod(t,24*3600*num_days) / 900)+1)]
	pv_infeed = t -> dem_data_twoweeks[!,:pv_infeed][Int(floor(mod(t,24*3600*num_days) / 900)+1)]

	demand_real_loads = t -> vcat(load(t), load(t), load(t), load(t)) ./ 20.
	demand_real = t -> vcat(load(t), pv_infeed(t)) ./ 50.
end

## Synthetic demand
begin
	struct demand_amp_var
		demand
	end

	function (dav::demand_amp_var)(t)
		index = Int(floor(t / (24*3600)))
		dav.demand[index + 1,:]
	end

    ## Generating synthetic demand function:
	demand_amp1 = demand_amp_var(repeat([80 80 80 10 10 10 40 40 40 40 40], outer=Int(N/4))')
	demand_amp2 = demand_amp_var(repeat([10 10 10 80 80 80 40 40 40 40 40], outer=Int(N/4))')
	demand_amp3 = demand_amp_var(repeat([60 60 60 60 10 10 10 40 40 40 40], outer=Int(N/4))')
	demand_amp4 = demand_amp_var(repeat([30 30 30 30 10 10 10 80 80 80 80], outer=Int(N/4))')
	demand_amp = t -> vcat(demand_amp1(t), demand_amp2(t), demand_amp3(t), demand_amp4(t))

	periodic_demand =  t-> demand_amp(t)./100 .* sin(t*pi/(24*3600))^2
	samples = 24*4
	inter = interpolate([.2 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
	residual_demand = t -> inter(1. + t / (24*3600) * samples)
	demand_synth = t -> periodic_demand(t) + residual_demand(t)
end
