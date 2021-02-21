# TODO: @inline

## Line constructor
function constructor(f::PowerLine)
	return StaticEdge(f! = f, dim = 1)
end

## Node constructors
function constructor(f::PV)
	# @assert f.η_gen(t) >= 0. && f.η_gen(t) <= 1. "Generation efficiency must be between 0 and 1"
	@assert f.M >= 0 "Inertia constant cannot be negative"
	@assert f.max_inverter >= 0 "Maximum inverter efficiency cannot be negative"
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_abs_LI])
end

function constructor(f::Load)
	# @assert f.η_gen(t) > 0. && f.η_gen(t) <= 1. "Load efficiency must be between 0 and 1"
	@assert f.M > 0 "Inertia constant cannot be negative"
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_abs_LI])
end
function constructor(f::Slack)
	# @assert f.η_gen(t) >= 0. && f.η_gen(t) <= 1. "Generation efficiency must be between 0 and 1"
	@assert f.M >= 0 "Inertia constant cannot be negative"
	@assert f.max >= 0 "Rated capacity cannot be negative"
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_abs_LI])
end

function constructor(f::Battery)
	# @assert f.η.gen(t) >= 0. && f.η.load(t) <= 1. "Generation and load efficiency must be between 0 and 1"
	@assert f.M >= 0 "Inertia constant cannot be negative"
	@assert f.max_inverter >= 0 "Maximum battery inverter output cannot be negative"
	@assert f.σ >= 0 "Battery self-discharge rate cannot be negative"
	return ODEVertex(f! = f, dim = 6, sym = [:ϕ, :ω, :χ, :integrated_LI, :level, :integrated_abs_LI])
end

function constructor(f::ThermalStorage)
	# @assert f.η.gen(t) >= 0. && f.η.load(t) <= 1. "Generation and load efficiency must be between 0 and 1"
	@assert f.M >= 0 "Inertia constant cannot be negative"
	@assert f.a_loss "Loss coefficient cannot be negative"
	return ODEVertex(f! = f, dim = 6, sym = [:ϕ, :ω, :χ, :integrated_LI, :level, :integrated_abs_LI])
end

function constructor(f::GenericNode)
	# @assert f.η.gen(t) >= 0. && f.η.load(t) <= 1. "Generation and load efficiency must be between 0 and 1"
	@assert f.M >= 0 "Inertia constant cannot be negative"
	@assert f.max_inverter >= 0 "Maximum battery inverter output cannot be negative"
	@assert f.σ >= "Battery self-discharge rate cannot be negative"
	return ODEVertex(f! = f, dim = 6, sym = [:ϕ, :ω, :χ, :integrated_LI, :level, :integrated_abs_LI])
end

# Create graph
@inline function Grid(n, l)
	N, L = length(n), length(l)
	nodes = [constructor(i) for i in n]
	lines = [constructor(i) for i in l]
    g = SimpleGraph(N)

	for i in 1:L
		add_edge!(g, l[i].from, l[i].to)
	end
	network_dynamics(nodes, lines, g)
end

# TODO: optional callback
# TODO: define types of attributes
# Simulate
@inline function simulate(nd, tspan, i0, solver, callback = nothing)
	cb_default = CallbackSet(PeriodicCallback(HourlyUpdate, l_hour),
					 		 PeriodicCallback(DailyUpdate, l_day))

	if callback !== nothing
		cb = CallbackSet(PeriodicCallback(HourlyUpdate, l_hour),
						 PeriodicCallback(DailyUpdate, l_day),
					     callback)
	else
		cb = cb_default
	end
	ode_prob = ODEProblem(nd, i0, tspan, callback = cb)
	sol = solve(ode_prob, solver)
	return sol
end
