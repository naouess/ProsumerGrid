"""
This code is a reproduction of the example impelemented in the paper
"Iterative learning control in prosumer-based microgrids with hierarchical control"
by Strenge et al. (2019)

# System of 4 nodes, graph degree 3
# Kappa = 1/h
# Coupling = 6
# Synthetic demand function (imported from "./demand/SyntheticDemand.jl/)
# 10 days

Plots are saved in the folder "./plots/Reprod_ILCpaper/"

"""
begin
	using NetworkDynamics, LightGraphs, Parameters, DifferentialEquations, GraphPlot
	using Interpolations
	using ToeplitzMatrices
	using DSP
	using Plots
	using LinearAlgebra
	using Statistics
	using Random
	Random.seed!(42)
end

# General parameters
begin
	num_days = 10
	l_day = 24*3600
	l_hour = 3600
	N = 4
end

begin
	dir = @__DIR__
	include("$dir/src/PowerNodes.jl")
	include("$dir/src/PowerLine.jl")
	include("$dir/src/Structs.jl")
	include("$dir/src/UpdateFunctions.jl")
	include("$dir/demand/SytheticDemand.jl")
end

# Parameters for ILC controller
begin
	# Number of ILC-Node
	vc = 1:N
	# Nodes with communication: here none
	cover = Dict([v => [] for v in vc])
	u = [zeros(1000,1);1;zeros(1000,1)];
	fc = 1/6
	a = digitalfilter(Lowpass(fc), Butterworth(2))
	Q1 = filtfilt(a, u) # Markov Parameter
	Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);
end

# Defining the nodes
dd = t->((periodic_demand(t) .+ residual_demand(t)))
begin
	myPV = PV(ξ = t -> dd(t)[1], η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T_inv = 1/0.04), M_inv = 1/5.)
	myLoad = Load(ξ = t -> dd(t)[2], η_load = t -> 1., LI = LI(kp = 110., ki = 0.004, T_inv = 1/0.045), M_inv = 1/4.8)
	mySlack = Slack(ξ = t -> dd(t)[3], η_gen = t -> 1., LI = LI(kp = 100., ki = 0.05, T_inv = 1/0.047), M_inv = 1/4.1)
	myLoad2 = Load2(ξ = t -> dd(t)[4], η_load = t -> 1., LI = LI(kp = 200., ki = 0.001, T_inv = 1/0.043), M_inv = 1/4.8)
	v1 = construct_vertex(myPV)
	v2 = construct_vertex(myLoad)
	v3 = construct_vertex(mySlack)
	v4 = construct_vertex(myLoad2)
end

begin
	# from, to parameters are not needed (or only for documentation and overview)
	myLine = PowerLine(from=1, to=2, K=6.)
	line = construct_edge(myLine)
end

begin
	nodes = [v1, v2, v3, v4]
	lines = [line for i in 1:6]
	g1 = random_regular_graph(N, 3)
	nd = network_dynamics(nodes, lines, g1)
end

# Defining ILC parameters
begin
	ILC_pars0= ILC(kappa = 1. /l_hour, mismatch_yesterday = zeros(24, N), daily_background_power = zeros(24, N),
	current_background_power = zeros(N), ilc_nodes = vc, ilc_cover = cover, Q = Q)
end

begin
	# Defining time span and initial values
	tspan = (0., num_days*l_day)
	ic = zeros(16)

	# Defining the callback function
	cb = CallbackSet(PeriodicCallback(HourlyUpdate(), l_hour), PeriodicCallback(DailyUpdate, l_day))

	# Passing first tuple element as parameters for nodes and nothing for lines
	ILC_p = (ILC_pars0, nothing)

	# Defining the ODE-Problem
	ode_problem = ODEProblem(nd, ic, tspan, ILC_p, callback = cb)
end

@time sol = solve(ode_problem, Rodas4())

## To plot state of system:
plot(sol, vars = syms_containing(nd, "ω"), legend = true)
plot(sol, vars = syms_containing(nd, "ϕ"), legend = true)
plot(sol, vars = syms_containing(nd, "integrated_LI"), legend = true)

## Extracting values out of solution found by solver
begin
	hourly_energy = zeros(24 * num_days + 1, N)
	for j = 1:N
		for i = 1:24*num_days+1
			hourly_energy[i, j] = sol((i-1)*3600)[4*j]
		end
	end

	# ILC power needed for a certain number of days
	# Also as norm and mean energy over each day

	# Values of ILC_power for day 0:
	ILC_power = zeros(num_days+2, 24, N)
	norm_energy_d = zeros(num_days,N)
	mean_energy_d = zeros(num_days,N)

	# Values for ILC_power for day 1:
	kappa = 1. /l_hour
	for j = 1:N
		ILC_power[2,:,j] = Q*(zeros(24,1) +  kappa * hourly_energy[1:24, j])
		norm_energy_d[1, j] = norm(hourly_energy[1:24, j])
		mean_energy_d[1, j] = mean(hourly_energy[1:24, j])
	end

	# Values for ILC_power from day 2 until end:
	for i= 2:num_days
		for j = 1:N
			ILC_power[i+1,:,j] .= Q*(ILC_power[i,:,j] .+  kappa*hourly_energy[(i-1)*24+1:i*24,j])
			norm_energy_d[i,j] = norm(hourly_energy[(i-1)*24+1:i*24,j])
			mean_energy_d[i,j] = mean(hourly_energy[(i-1)*24+1:i*24,j])
		end
	end
end

# plot(hourly_energy[:, 2] ./ maximum(hourly_energy[:, 2]))
# plot(ILC_power[10,:,1])

## Plots
using LaTeXStrings

# NODE WISE second-wisenode = 1
begin
	node = 1
	p1 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[node], alpha=0.2, label = latexstring("P^d_$node"),linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600,hourly_energy[1:num_days*24, node]./3600, label=latexstring("y_$node^{c,h}"),linewidth=3)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_node[1:num_days*24], label=latexstring("\$u_$node^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	      xtickfontsize=14, legendfontsize=10, linewidth=3, yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
	title!(L"j = 1")
	savefig("$dir/plots/Reprod_ILCpaper/demand_seconds_node_$(node)_hetero.png")
end

begin
	node = 2
	p2 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[node], alpha=0.2, label = latexstring("P^d_$node"),linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600,hourly_energy[1:num_days*24,node]./3600, label=latexstring("y_$node^{c,h}"),linewidth=3)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_node[1:num_days*24], label=latexstring("\$u_$node^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
				   xtickfontsize=14, yticks=false, legendfontsize=10, linewidth=3,legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
	title!(L"j = 2")
	savefig("$dir/plots/Reprod_ILCpaper/demand_seconds_node_$(node)_hetero.png")
end

begin
	node = 3
	p3 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[node], alpha=0.2, label = latexstring("P^d_$node"),linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600,hourly_energy[1:num_days*24,node]./3600, label=latexstring("y_$node^{c,h}"),linewidth=3)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_node[1:num_days*24], label=latexstring("\$u_$node^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	               xtickfontsize=18, legendfontsize=10, linewidth=3,xaxis=("days [c]",font(14)),yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
	title!(L"j = 3")
	savefig("$dir/plots/Reprod_ILCpaper/demand_seconds_node_$(node)_hetero.png")
end

begin
	node = 4
	p4 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[node], alpha=0.2, label = latexstring("P^d_$node"),linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600,hourly_energy[1:num_days*24,node]./3600, label=latexstring("y_$node^{c,h}"),linewidth=3)#, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_node[1:num_days*24], label=latexstring("\$u_$node^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	               xtickfontsize=18,
	    		   legendfontsize=10, yticks=false, linewidth=3,xaxis=("days [c]",font(14)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-0.7,1.5)
	title!(L"j = 4")
	savefig("$dir/plots/Reprod_ILCpaper/demand_seconds_node_$(node)_hetero.png")
end

l = @layout [a b; c d]
plot_demand = plot(p1,p2,p3,p4,layout = l)
savefig(plot_demand, "$dir/plots/Reprod_ILCpaper/demand_seconds_all_nodes_hetero.png")

psum = plot()
ILC_power_hourly_mean_sum = vcat(ILC_power[:,:,1]'...) .+ vcat(ILC_power[:,:,2]'...) .+ vcat(ILC_power[:,:,3]'...) .+ vcat(ILC_power[:,:,4]'...)
plot!(0:num_days*l_day, t -> (dd(t)[1] .+ dd(t)[2] .+ dd(t)[3] .+ dd(t)[4]), alpha=0.2, label = latexstring("\$P^d_j\$"),linewidth=3, linestyle=:dot)
plot!(1:3600:24*num_days*3600,(hourly_energy[1:num_days*24,1] + hourly_energy[1:num_days*24,2] + hourly_energy[1:num_days*24,3] + hourly_energy[1:num_days*24,4])./3600, label=latexstring("y_j^{c,h}"),linewidth=3, linestyle=:dash)
plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], label=latexstring("\$u_j^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
               xtickfontsize=18,legend=false, legendfontsize=10, linewidth=3,xaxis=("days [c]",font(14)),  yaxis=("normed power",font(14)),lc =:black, margin=5Plots.mm)
# ylims!(-0.7,1.5)
# title!("Initial convergence")
savefig(psum,"$dir/plots/Reprod_ILCpaper/demand_seconds_sum_hetero.png")
