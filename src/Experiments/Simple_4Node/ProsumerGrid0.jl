begin
	using NetworkDynamics, LightGraphs, Parameters, DifferentialEquations, GraphPlot
	using Interpolations
	using Sundials # to use the solver CVODE_BDF()
	using ODEInterfaceDiffEq # to use the solver radau()
	using ToeplitzMatrices
	using DSP
	using Plots
	using LinearAlgebra
	using Statistics
	using Random
	Random.seed!(42)
end

#=
Minimal example to test out the defined nodes and Powerline
=#

# General parameters
begin
	num_days = 10
	l_day = 24*3600
	l_hour = 3600
	N = 4
end

begin
	dir = @__DIR__
	include("$dir/PowerNodes0.jl")
	include("$dir/PowerLine0.jl")
	include("$dir/Structs.jl")
	include("$dir/UpdateFunctions.jl")
	include("$dir/Demand.jl")
end

# Parameters for ILC controller
begin
	# Number of ILC-Node
	vc = 1:N
	# Nodes with communication: here none
	cover = Dict([v => [] for v in vc]) # cover = Dict([1 => []])
	u = [zeros(1000,1);1;zeros(1000,1)];
	fc = 1/6
	a = digitalfilter(Lowpass(fc), Butterworth(2))
	Q1 = filtfilt(a, u) # Markov Parameter
	Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);
end

begin
	nodes = [node, node, node, node]
	lines = [line, line, line, line, line, line]
	g1 = random_regular_graph(N, 3)
	# gplot(g1, nodelabel=1:4)
	nd = network_dynamics(nodes, lines, g1)
end

# Defining ILC parameters
begin
	ILC_pars1 = ILC(kappa = 1/l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
end

# Defining leaky integrator parameters
begin
	LI_pars_PV = LI(kp = 400., ki = 0.05, T_inv = 1/0.04)
	LI_pars_load = LI(kp = 110., ki = 0.004, T_inv = 1/0.045)
	LI_pars_Slack = LI(kp = 100., ki = 0.05, T_inv = 1/0.047)
	LI_pars_Load2 = LI(kp = 200., ki = 0.001, T_inv = 1/0.043)
end

# Defining compound parameters
begin
	par_PV = compound_pars(ξ = t -> (periodic_demand(t)+residual_demand(t))[1], ILC=ILC_pars1, LI=LI_pars_PV, M_inv = 1/5.)
	par_load = compound_pars(ξ = t -> (periodic_demand(t)+residual_demand(t))[2], ILC=ILC_pars1, LI=LI_pars_load, M_inv = 1/4.8)
	par_slack = compound_pars(ξ = t -> (periodic_demand(t)+residual_demand(t))[3], ILC=ILC_pars1, LI=LI_pars_Slack, M_inv = 1/4.1)
	par_load2 = compound_pars(ξ = t -> (periodic_demand(t)+residual_demand(t))[4], ILC=ILC_pars1, LI=LI_pars_Load2, M_inv = 1/4.8)
	# par_PV = compound_pars(ξ = t -> 2., ILC=ILC_pars1, LI=LI_pars_PV, M_inv = 1/5.)
	# par_load = compound_pars(ξ = t -> -0.8, ILC=ILC_pars1, LI=LI_pars_load, M_inv = 1/4.8)
	# par_slack = compound_pars(ξ = t -> 1., ILC=ILC_pars1, LI=LI_pars_Slack, M_inv = 1/4.1)
	# par_load2 = compound_pars(ξ = t -> -1., ILC=ILC_pars1, LI=LI_pars_Load2, M_inv = 1/4.8)
end

plot(1:24*3600*7, t -> par_PV.ξ(t))
plot(1:24*3600*7, t -> par_load.ξ(t))
# plot(1:24*3600*7, t -> par_load2.ξ(t))
# plot(1:24*3600*7, t -> par_slack.ξ(t))

begin
	# Defining time span and initial values
	tspan = (0., num_days*l_day)
	ic = ones(16)

	# Defining the callback function
	cb = CallbackSet(PeriodicCallback(HourlyUpdate(), l_hour), PeriodicCallback(DailyUpdate, l_day))

	# Passing first tuple element as parameters for nodes and nothing for lines
	parameters = ([par_PV, par_load, par_slack, par_load2], nothing) #

	# Defining the ODE-Problem
	ode_problem = ODEProblem(nd, ic, tspan, parameters, callback = cb)
end

@time sol = solve(ode_problem, Rodas4())
# Rosenbrock23()

## Extracting values out of solution found by solver

hourly_energy = zeros(24 * num_days + 1, N)
hourly_theta = zeros(24 * num_days + 1, N)
for j = 1:N
	for i = 1:24*num_days+1
		hourly_energy[i, j] = sol((i-1)*3600)[4*j]
		# hourly_theta[i, j] = sol((i-1)*3600)[4*j - 2]
	end
end
plot(hourly_energy)
plot(hourly_theta)

begin
	# ILC power needed for a certain number of days
	# Also as norm and mean energy over each day
	ILC_power = zeros(num_days+2, 24, N)
	norm_energy_d = zeros(num_days,N)
	mean_energy_d = zeros(num_days,N)

	kappa = 1/l_hour
	for j = 1:N
		ILC_power[2,:,j] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,j])
		norm_energy_d[1,j] = norm(hourly_energy[1:24,j])
		mean_energy_d[1,j] = mean(hourly_energy[1:24,j])
	end

	for j = 1:N
		for i= 2:num_days
			ILC_power[i+1,:,j] .= Q*(ILC_power[i,:,j] .+  kappa*hourly_energy[(i-1)*24+1:i*24,j])
			norm_energy_d[i,j] = norm(hourly_energy[(i-1)*24+1:i*24,j])
			mean_energy_d[i,j] = mean(hourly_energy[(i-1)*24+1:i*24,j])
		end
	end

	ILC_power_agg = maximum(mean(ILC_power.^2,dims=3),dims=2)
	ILC_power_agg = mean(ILC_power,dims=2)
	ILC_power_agg_norm = norm(ILC_power)
	ILC_power_hourly = vcat(ILC_power[:,:,1]'...)
end



## Plots
using LaTeXStrings

begin
# hourly plotting
	plot(1:num_days*24, ILC_power_hourly[1:24*num_days] , legend=:topleft, label=L"$ u_j^{ILC}$", ytickfontsize=14,
	               xtickfontsize=18,
	    		   legendfontsize=12, linewidth=3,xaxis=("time [h]",font(14)), yaxis=("normed power",font(14)))
	plot!(1:num_days*24+1,mean(hourly_energy, dims=2)/3600 , label=L"y^{c,h}", linewidth=3)
end

# second-wise
begin
	plot(1:3600:num_days*24*3600,  ILC_power_hourly[1:num_days*24]./ maximum(ILC_power_hourly), label=L"$P_{ILC, j}$", ytickfontsize=14,
	               xtickfontsize=18,
	    		   legendfontsize=10, linewidth=3,xaxis=("time [s]",font(14)), yaxis=("normed quantities [a.u.]",font(14)))
	plot!(1:3600:24*num_days*3600,mean(hourly_energy[1:num_days*24], dims=2) ./ maximum(hourly_energy), label=L"y_h",linewidth=3, linestyle=:dash)
end
