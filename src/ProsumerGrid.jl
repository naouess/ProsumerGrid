begin
	using NetworkDynamics, LightGraphs, Parameters, DifferentialEquations, GraphPlot
	using Sundials # to use the solver CVODE_BDF()
	using ODEInterfaceDiffEq # to use the solver radau()
	using ToeplitzMatrices
	using DSP
	using Plots
end

begin
	dir = @__DIR__
	include("$dir/PowerNodes.jl")
	include("$dir/PowerLine.jl")
	include("$dir/structs.jl")
	include("$dir/UpdateFunctions.jl")
end

#=
Minimal example to test out the defined nodes and Powerline
=#

# General parameters
begin
	num_days = 35
	l_day = 24*3600
	l_hour = 3600
	N = 4
end

# Parameters for ILC controller
begin
	# Number of ILC-Node (here: without communication)
	vc = 1:N
	cover = Dict([v => [] for v in vc])
	u = [zeros(1000,1);1;zeros(1000,1)];
	fc = 1/6
	a = digitalfilter(Lowpass(fc), Butterworth(2))
	Q1 = filtfilt(a, u) # Markov Parameter
	Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);
end

# Defining the nodes
begin
	myPV = PV(ξ = t -> 1000., η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T = 0.04), M = 5.)
	myLoad = Load(ξ = t -> -1000., η_load = t -> 1., LI = LI(kp = 200., ki = 0.004, T = 0.045), M = 4.8)
	mySlack = Slack(ξ = t -> 800., η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T = 0.04), M = 4.1)
	myLoad2 = Load(ξ = t -> -800., η_load = t -> 1., LI = LI(kp = 200., ki = 0.004, T = 0.045), M = 4.8)
	v1 = construct_vertex(myPV)
	v2 = construct_vertex(myLoad)
	v3 = construct_vertex(mySlack)
	v4 = construct_vertex(myLoad2)
	nodes = [v1, v2, v3, v4]
end

begin
	# from, to parameters are not needed (or only for documentation and overview)
	myLine = PowerLine(from=1, to=2, K=6)
	l1 = construct_edge(myLine)
	lines = [l1, l1, l1, l1]
end

# begin
#     g1 = random_regular_graph(length(nodes), 3)
#     gplot(g1, nodelabel=1:4)
# end

begin
	g = SimpleGraph(4)
	add_edge!(g, 1, 2)
	add_edge!(g, 1, 4)
	add_edge!(g, 3, 4)
	add_edge!(g, 3, 2)
	gplot(g, nodelabel=1:4)
end

nd = network_dynamics(nodes, lines, g)

tspan = (0., num_days * l_day)
ic = ones(5 * N)


ILC_pars1 = ILC(kappa = 1/l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)

# Defining the callback function
cb = CallbackSet(PeriodicCallback(HourlyUpdate(), l_hour), PeriodicCallback(DailyUpdate_X, l_day))
# Passing first tuple element as parameters for nodes and nothing for lines
ILC_p = ([ILC_pars1, ILC_pars1, ILC_pars1, ILC_pars1], nothing)
# Defining the ODE-Problem
ode_problem = ODEProblem(nd, ic, tspan, ILC_p, callback = cb)

@time sol = solve(ode_problem, Rodas4())
# Other solvers to try out: Rosenbrock23(), Rodas4(autodiff=false)
# CVODE_BDF() doesn't use mass matrices
# radau() doesn't have DEOptions (?)

####
# Extracting values out of solution found by solver
hourly_energy = zeros(24 * num_days + 1, N)
for i = 1:24*num_days+1
	for j = 1:N
		hourly_energy[i, j] = sol((i-1)*3600)[4*j]
	end
end
plot(hourly_energy)

# ILC power needed for a certain number of days
ILC_power = zeros(num_days+2, 24, N)

kappa = 1/l_hour
for j = 1:N
	ILC_power[2,:,j] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,j])
end

for i= 2:num_days
	for j = 1:N
		ILC_power[i+1,:,j] = Q*(ILC_power[i,:,j] +  kappa*hourly_energy[(i-1)*24+1:i*24,j])
	end
end
