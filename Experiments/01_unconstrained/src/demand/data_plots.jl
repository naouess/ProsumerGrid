using JLD2, FileIO, GraphIO, CSV, DataFrames
using Statistics
using StatsBase
using Plots
using Interpolations

num_days = 10
l_day = 3600*24
N = 4

begin
	## Generating demand function out of real minigrid data
	dem_data_twoweeks = CSV.read("$dir/src/demand/TwoWeeks.csv")

	load = t->dem_data_twoweeks[!,:load][Int(floor(mod(t,24*3600*num_days) / 900)+1)]
	pv_infeed = t->dem_data_twoweeks[!,:pv_infeed][Int(floor(mod(t,24*3600*num_days) / 900)+1)]

	demand_real = t -> vcat(load(t), pv_infeed(t)) ./ 50.
end

# # Plotting
using Random
Random.seed!(42)

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

begin
	# plot exact power balance
	p_demandreal= plot()
	plot!(1:24*num_days*3600, t-> demand_real(t)[1] .* (50), linecolor=:lightblue,
	fill=(0, 0, :lightblue),
	xticks=(0:3600:3600*24*num_days, string.(0:24*num_days)),
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=6, legendfontsize=10,
	xaxis=("Hours",font(10)),  yaxis=("Energy demand in kWh",font(10)), margin=5Plots.mm,
	)
	# savefig("$dir/plots/real_demand_beautiful_oneday.pdf")
end

begin
	# plot exact power balance
	p_demandreal= plot()
	plot!(1:24*num_days*3600, t-> demand_real(t)[2], linecolor=:lightblue,
	fill=(0, 0, :lightblue),
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	# alpha=1,
	#linestyle=:dot,
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm,
	# ribbon=mismatch_solt
	)
	savefig("$dir/plots/real_infeed_beautiful.pdf")
end

begin
	# plot exact power balance
	p_demandreal= plot()
	plot!(1:24*num_days*3600, t-> demand_real(t)[1] .- demand_real(t)[2], linecolor=:lightblue,
	fill=(0, 0, :lightblue),
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	# alpha=1,
	#linestyle=:dot,
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm,
	# ribbon=mismatch_solt
	)
	savefig("$dir/plots/real_net_beautiful.pdf")
end

data_load = [demand_real(i)[1] for i in 0:3600:24*num_days*3600-1]
data_pv = [demand_real(i)[2] for i in 0:3600:24*num_days*3600-1]
data_net = [demand_real(i)[1] .- demand_real(i)[2] for i in 0:3600:24*num_days*3600-1]
a = autocor(data_load, demean = true,[i for i in 0:num_days*24-1])
b = autocor(data_pv, demean = true , [i for i in 0:num_days*24-1])
d = autocor(data_net, demean = true , [i for i in 0:num_days*24-1])
begin
	p2 = plot(bar(a), xticks = (0:24:num_days*24, string.(0:num_days)), legend=false, fillcolor=:turquoise, linecolor=:grey,
		 xaxis=("Days", font(10)), yaxis=("ACF", font(10)), margin=5Plots.mm)
	ylims!(-0.6, 1.)
	# savefig("$dir/plots/autocor_load.pdf")
end

begin
	p1 = plot(bar(b), xticks = (0:24:num_days*24, string.(0:num_days)), legend=false, fillcolor=:turquoise, linecolor=:grey,
		 xaxis=("Days", font(10)), yaxis=("ACF", font(10)), margin=5Plots.mm)
	ylims!(-0.6, 1)
	# savefig("$dir/plots/autocor_infeed.pdf")
end

begin
	p4 = plot(bar(d), xticks = (0:24:num_days*24, string.(0:num_days)), legend=false, fillcolor=:turquoise, linecolor=:grey,
		 xaxis=("Days", font(10)), yaxis=("ACF", font(10)), margin=5Plots.mm)
	ylims!(-0.6, 1)
	savefig("$dir/plots/autocor_net_real.pdf")
end

data_synth = [(demand_synth(i)[1] + demand_synth(i)[2] + demand_synth(i)[3] + demand_synth(i)[4]) for i in 0:3600:24*num_days*3600-1]
data_synth_oneload = [demand_synth(i)[1] for i in 0:3600:24*num_days*3600-1]
data_periodic = [sum(periodic_demand(i)) for i in 0:3600:24*num_days*3600-1]
c = autocor(data_synth, demean = true , [i for i in 0:num_days*24-1])
c_1 = autocor(data_periodic, demean = true , [i for i in 0:num_days*24-1])
f = autocor(data_synth_oneload, demean = true , [i for i in 0:num_days*24-1])
begin
	p3= plot(bar(c), xticks = (0:24:num_days*24, string.(0:num_days)), legend=false, fillcolor=:turquoise, linecolor=:grey,
		 xaxis=("Days", font(10)), yaxis=("ACF", font(10)), margin=5Plots.mm)
	ylims!(-1, 1.)
	savefig("$dir/plots/autocor_synth_aggregated.pdf")
end

begin
	p5= plot(bar(f), xticks = (0:24:num_days*24, string.(0:num_days)), legend=false, fillcolor=:turquoise, linecolor=:grey,
		 xaxis=("Days", font(10)), yaxis=("ACF", font(10)), margin=5Plots.mm)
	ylims!(-1, 1.)
	savefig("$dir/plots/autocor_synth_oneload.pdf")
end

begin
	p5= plot(bar(c_1), xticks = (0:24:num_days*24, string.(0:num_days)), legend=false, fillcolor=:turquoise, linecolor=:grey,
		 xaxis=("Days", font(10)), yaxis=("ACF", font(10)), margin=5Plots.mm)
	ylims!(-1, 1.)
	savefig("$dir/plots/autocor_synth_periodic.pdf")
end
