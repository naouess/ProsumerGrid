begin
	using NetworkDynamics, LightGraphs, Parameters, DifferentialEquations, GraphPlot
	using Interpolations
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
	include("$dir/Structs.jl")
	include("$dir/UpdateFunctions.jl")
	# include("$dir/Demand.jl")
end

#=
Minimal example to test out the defined nodes and Powerline
=#

# General parameters
begin
	num_days = 20
	l_day = 24*3600
	l_hour = 3600
	N = 2
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

# Defining the nodes
begin
	myPV = PV(ξ = t -> 2., η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T_inv = 1/0.04), M_inv = 1/5.) #(periodic_demand(t)+residual_demand(t))/N
	myLoad = Load(ξ = t -> -1.8, η_load = t -> 1., LI = LI(kp = 110., ki = 0.004, T_inv = 1/0.045), M_inv = 1/4.8)
	# mySlack = Slack(ξ = t -> 0.8, η_gen = t -> 1., LI = LI(kp = 100., ki = 0.05, T_inv = 1/0.047), M_inv = 1/4.1)
	# myLoad2 = Load(ξ = t -> -0.8, η_load = t -> 1., LI = LI(kp = 200., ki = 0.001, T_inv = 1/0.043), M_inv = 1/4.8)
	v1 = construct_vertex(myPV)
	v2 = construct_vertex(myLoad)
	# v3 = construct_vertex(mySlack)
	# v4 = construct_vertex(myLoad2)
	# nodes = [v1, v2, v3, v4]
	nodes = [v1, v2]
end

begin
	# from, to parameters are not needed (or only for documentation and overview)
	myLine = PowerLine(from=1, to=2, K=6)
	l1 = construct_edge(myLine)
	lines = [l1]
	# lines = [l1, l1, l1, l1, l1, l1]
end

# begin
#     g1 = random_regular_graph(length(nodes), 1)
#     gplot(g1, nodelabel=1:2)
# end

begin
	g = SimpleDiGraph(2)
	add_edge!(g, 1, 2)
	# add_edge!(g, 1, 4)
	# add_edge!(g, 3, 4)
	# add_edge!(g, 3, 2)
	gplot(g, nodelabel=1:2)
end

nd = network_dynamics(nodes, lines, g)

# Defining time span and initial values
tspan = (0., num_days*l_day)
# ic = [1.175, 0., 1., 1., 0., 1., 0., -1., -1., 0., 13.7/12, 0., 0.8, 0.8, 0., 6.1/6, 0., -0.8, -0.8, 0.] # zeros(5 * N)
ic = [7/6, 0., 1., 1., 0., 1., 0., -1., -1., 0.]
# ic = zeros(10)

# Defining ILC parameters
ILC_pars1 = ILC(kappa = 1/l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)

# Defining the callback function
cb = CallbackSet(PeriodicCallback(HourlyUpdate(), l_hour), PeriodicCallback(DailyUpdate, l_day))

# Passing first tuple element as parameters for nodes and nothing for lines
# ILC_p = ([ILC_pars1, ILC_pars1, ILC_pars1, ILC_pars1], nothing)
ILC_p = ([ILC_pars1, ILC_pars1], nothing)

# Defining the ODE-Problem
ode_problem = ODEProblem(nd, ic, tspan, ILC_p, callback = cb)

@time sol =  solve(ode_problem, Rodas4()) # @debug
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
