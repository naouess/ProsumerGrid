### Code investigating possibilities of haing lower-layer conrol on specific nodes only ####
begin
	dir = @__DIR__
	include("$dir/src/ProsumerGrid.jl")
end

begin
	using ToeplitzMatrices
	using DSP
	using LinearAlgebra
end

# General parameters
begin
	num_days = 10
	l_day = 24*3600
	l_hour = 3600
	N = 4
end

# Import demand functions
include("$dir/src/MinigridDemand.jl")
d = t -> demand_real(t)
# Define ILC parameters for each node
begin
	vc = 1:N
	cover = Dict([v => [] for v in vc])

	# No Q-filter applied
	Q = Matrix(1.0I, 24, 24)

	ILC_pars1 = ILC(kappa = 1/l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
				    current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars2 = ILC(kappa = 1 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
					current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars3 = ILC(kappa = 1 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
					current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars4 = ILC(kappa = 1 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
					current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
end

# Define LI parameters for each node
begin
	# LI_pars1 = LI(kp = 400., ki = 0.05, T = 0.04)
	# No LI at load node:
	LI_pars1 = LI(kp = 100., ki = 0., T = 0.)
	# LI_pars2 = LI( kp = 110., ki = 0.004, T = 0.045)
	# No LI at load node:
	LI_pars2 = LI( kp = 100., ki = 0., T = 0.)
	LI_pars3 = LI(kp = 100., ki = 0.05, T = 0.047)
	LI_pars4 = LI(kp = 200., ki = 0.001, T = 0.043)
end

# Define nodes
begin
	Node1 = PV(ξ = t -> d(t)[2], η_gen = t-> 1., M = 5, LI = LI_pars1, ILC = ILC_pars1)
	Node2 = Load(ξ = t -> d(t)[1], η_load = t-> 1., M = 4.8, LI = LI_pars2, ILC = ILC_pars2)
	Node3 = Slack(η_gen = t -> 1., M = 4.8, LI = LI_pars3, ILC = ILC_pars3)
	Node4 = Battery(η = η(gen = t -> 1., load=t->1.), M = 4.8, C= 310 / 50, LI = LI_pars4, ILC = ILC_pars4)

	nodes = [Node1, Node2, Node3, Node4]
end

# Define lines
begin
	Line1 = PowerLine(from=1, to=2, K=6.)
	Line2 = PowerLine(from=1, to=4, K=6.)
	Line3 = PowerLine(from=3, to=2, K=6.)
	Line4 = PowerLine(from=4, to=2, K=6.)
	Line5 = PowerLine(from=3, to=4, K=6.)

	lines = [Line1, Line2, Line3, Line4, Line5]
end

# Construct the grid
pg = Grid(nodes, lines)

# Plot the created graph
# using GraphPlot
# gplot(pg.f.graph, nodelabel = 1:length(nodes))

# Define the time span
tspan = (0., num_days * 24 * 3600) # 10 days of simulation

# Define initial values
i0 = zeros(21)

# Solving the problem
sol = simulate(pg, tspan, i0, Rodas4())

# Quick plots
using Plots
plot(sol, vars = syms_containing(pg, "ω"), legend = true)
plot(sol, vars = syms_containing(pg, "ϕ"), legend = true)
plot(sol, vars = syms_containing(pg, "integrated_LI"), legend = true)
plot(sol, vars = syms_containing(pg, "level"), legend = true)
