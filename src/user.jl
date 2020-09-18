## TODO
# Add dependencies
# using ProsumerGrid
include("ProsumerGrid.jl")

# TODO:
# (done) generalize p.current_background_power[]:
# (done) add w component
# (done) try with mass_matrix
# (done) add battery model
# (done) add efficiencies
# (    ) add bounds
# (    ) add generic node model for experimentations :)

begin
	using Interpolations
	using ToeplitzMatrices
	using DSP
	using Plots
	using LinearAlgebra
	using Random
	Random.seed!(42)
	using GraphPlot
	# using Sundials # to use the solver CVODE_BDF()
	# using ODEInterfaceDiffEq # to use the solver radau()
end

# General parameters
begin
	num_days = 7
	l_day = 24*3600
	l_hour = 3600
	N = 4
end

# Import demand functions
begin
	dir = @__DIR__
	include("$dir/MinigridDemand.jl")
end

# Define ILC parameters for each node
begin
	# ILC-Node
	vc = 1:N
	# Nodes with communication: here none
	cover = Dict([v => [] for v in vc])
	u = [zeros(1000,1);1;zeros(1000,1)];
	fc = 1/6
	a = digitalfilter(Lowpass(fc), Butterworth(2))
	Q1 = filtfilt(a, u)
	Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);

	ILC_pars0 = ILC(kappa = 1. /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars1 = ILC(kappa = 1. /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars2 = ILC(kappa = 1. /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars3 = ILC(kappa = 1. /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
end

# Define nodes, lines and power grid structure
begin
	dd = t -> demand(t)
	myPV = PV(ξ = t -> dd(t)[2], η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T_inv = 1/0.04), ILC = ILC_pars0, M_inv = 1/5.)
	myLoad = Load(ξ = t -> dd(t)[1], η_load = t -> 1., LI = LI(kp = 110., ki = 0.004, T_inv = 1/0.045), ILC = ILC_pars1, M_inv = 1/4.8)
	mySlack = Slack(ξ = t -> 0., η_gen = t -> 1., LI = LI(kp = 100., ki = 0.05, T_inv = 1/0.047), ILC = ILC_pars2, M_inv = 1/4.1)
	myBattery = Battery(η = η(gen=t->1., load=t->1.), LI = LI(kp = 200., ki = 0.001, T_inv = 1/0.043), ILC = ILC_pars3, M_inv = 1/4.8, C=310. / 30.) # C in kWh, normed
	# myLoad2 = Load(ξ = t -> abs.(dd(t)[4]), η_load = t -> 1., LI = LI(kp = 200., ki = 0.001, T_inv = 1/0.043), ILC = ILC_pars3, M_inv = 1/4.8)

	nodes = [constructor(myPV), constructor(myLoad), constructor(mySlack), constructor(myBattery)]

	myLine1 = PowerLine(from=1, to=2, K=6.)
	myLine2 = PowerLine(from=1, to=4, K=6.)
	myLine3 = PowerLine(from=3, to=2, K=6.)
	myLine4 = PowerLine(from=4, to=2, K=6.)
	# myLine5 = PowerLine(from=1, to=2, K=6.)
	mylines = [myLine1, myLine2, myLine3, myLine4] #, myLine5]
	lines = [StaticEdge(f! = myLine, dim = 1) for myLine in mylines]

	g = SimpleGraph(N)
	add_edge!(g, 1, 2)
	add_edge!(g, 1, 4)
	add_edge!(g, 3, 2)
	# add_edge!(g, 3, 4)
	add_edge!(g, 2, 4)

	nd = network_dynamics(nodes, lines, g)
end
gplot(g, nodelabel=1:4, edgelabel = 1:5)

begin
	tspan = (0., num_days*l_day)
	# i0 = [0., 0., 0., -dd(0)[1], 0., 0., 0., 0., dd(0)[2], 0., 0., 0., 0., 0., 0., 0., 0., dd(0)[4], 0.]
	i0 = ones(21) * 0.5

	# Define the SavingCallback:
	# saved_edgevalues = SavedValues(Float64, Array{Float64, 1})
	# save_callback = SavingCallback(save_edges, saved_edgevalues)

	# Define Periodic and Bounds watching Callbacks
	cb = CallbackSet(PeriodicCallback(HourlyUpdate, l_hour),
					 PeriodicCallback(DailyUpdate, l_day))#,
			 		 # watch_storage_1(nd))#,
					 # watch_storage_2(nd))
			 		 # watch_storage_discharge(nd), watch_storage_charge(nd)) # , initial_affect= false)

	# Define the ODE-Problem
	ode_problem = ODEProblem(nd, i0, tspan, callback = cb)#, reltol = 1e-3, abstol = 1e-3, dtmax = 0.01)
end

@time sol = solve(ode_problem, Rodas4())
@show sol.destats

# Quick plots for validation
plot(sol, vars = syms_containing(nd, "ω"), legend = true)
# ylims!(-0.002, 0.002)
plot(sol, vars = syms_containing(nd, "ϕ"), legend = true)
plot(sol, vars = syms_containing(nd, "integrated_LI"), legend = true)
# plot(sol, vars = syms_containing(nd, "curtailed"), legend = true)
plot(sol, vars = syms_containing(nd, "level"), legend = true)

#= More plots
plot_flow(saved_edgevalues, nd)
indices1 = idx_containing(nd, :curtailed)
vars = syms_containing(nd, :curtailed)

curtailed = zeros(24 * num_days + 1, length(vars))
for j = 1:length(vars)
	for i = 1:24*num_days+1
		curtailed[i, j] = sol((i-1)*3600)[indices[j]]
	end
end
plot(curtailed ./ 3600)
=#

# Extract hourly control powers
hourly_LI, ILC_power = hourly_energy(sol, nd, num_days, N)
plot(hourly_LI)

# Plot battery stuff
begin
	indices1 = idx_containing(nd, :level)
	vars = syms_containing(nd, :level)

	level = zeros(24 * num_days + 1, length(vars))
	for j = 1:length(vars)
		for i = 1:24*num_days+1
			level[i, j] = sol((i-1)*3600)[indices1[j]]
		end
	end
end

plot(
    [level[:, 1], hourly_LI[:, 4]],
    # label = ["height" "spring length"],
    layout = (2, 1),
)

# using LaTeXStrings

# Nodewise plots
begin
	node = 1
	p1 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[2], alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, hourly_LI[1:num_days*24, node]./3600,linewidth=3)
	plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	      xtickfontsize=14, legendfontsize=10, linewidth=3, yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
	# savefig("$dir/plots/demand_seconds_node_$(node)_hetero.png")
end

begin
	node = 2
	p2 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[1], alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, hourly_LI[1:num_days*24, node]./3600,linewidth=3)
	plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	      xtickfontsize=14, legendfontsize=10, linewidth=3, yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
    ylims!(-0.7,1.5)
end

begin
	node = 3
	p3 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> 0, alpha=0.2,linewidth=3, linestyle=:dot)
	# plot!(1:3600:num_days*l_day, curtailed[1:num_days*24]./3600, alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, hourly_LI[1:num_days*24, node]./3600,linewidth=3)
	plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	      xtickfontsize=14, legendfontsize=10, linewidth=3, yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
end
begin
	node = 4
	p4 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> 0., alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, hourly_LI[1:num_days*24, node]./3600, linewidth=3)
	plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	      xtickfontsize=14, legendfontsize=10, linewidth=3, yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
end

begin
	l = @layout [a b; c d]
	plot_demand = plot(p1,p2,p3,p4, layout = l)
end

begin
	psum = plot()
	ILC_power_hourly_mean_sum =    vcat(ILC_power[:,:,1]'...) .+ vcat(ILC_power[:,:,2]'...) .+ vcat(ILC_power[:,:,3]'...) .+ vcat(ILC_power[:,:,4]'...)
	plot!(0:num_days*l_day, t -> (dd(t)[1] .- dd(t)[2]), alpha=0.2,linewidth=3, linestyle=:dot)
	# plot!(1:3600:num_days*l_day, curtailed[1:num_days*24]./3600, alpha=0.2, label = latexstring("Curtailed_$node"),linewidth=3, linestyle=:dot, lc = :red)
	plot!(1:3600:24*num_days*3600, (hourly_LI[1:num_days*24,1] + hourly_LI[1:num_days*24,2] + hourly_LI[1:num_days*24,3] + hourly_LI[1:num_days*24,4])./3600,linewidth=3, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	               xtickfontsize=18,legend=false, legendfontsize=10, linewidth=3,xaxis=("days [c]",font(14)),  yaxis=("normed power",font(14)),lc =:black, margin=5Plots.mm)
end
