module ProsumerGrid

begin
	using Parameters
	using Pkg
	# ? Why use this?
	Pkg.instantiate()
	using PowerDynamics
	using NetworkDynamics

	# TODO remove unneeded imports
	import PowerDynamics: dimension, symbolof, construct_vertex, construct_edge
end

begin
	dir = @__DIR__
	include("$dir/nodes/NodeTypesPD.jl")
	include("$dir/lines/PowerLine.jl")
end

@with_kw mutable struct η
	gen
	load
	"""
	Constructor
	"""
	function η(gen, load)
			new(gen, load)
	end
end # TODO: Define outer constructors for η for the cases when only gen or load is defined

@with_kw mutable struct u
	ILC
	LI
	"""
	Constructor
	"""
	function u(ILC, LI)
			new(ILC, LI)
	end
end

#=
Minimal example to test out the defined nodes and Powerline
No ILC solving integrated yet. u is set as a parameter.
=#

u_PV = u(800., 0.)
u_load = u(-900.,0.)
u_slack = u(100.,0.)

myPV = PV(ξ=800., η_gen=1., LI = u_PV.LI, ILC = u_PV.ILC, M=10.)
myLoad = Load(ξ=-900., η_load = 1., LI = u_load.LI, ILC = u_load.ILC, M=10.)
#mySlack = Slack(ξ=100., η_gen=1., LI = u_slack.LI, ILC = u_slack.ILC, M=1.)
#slack2 = SlackAlgebraic(U=1.) # To avoid triggering warnings
#Battery = Battery()

nodes = [myPV, myLoad]#, mySlack] #, slack2]

lines=[PowerLine(from=1,to=2,K=100)]#, PowerLine(from=3,to=2,K=6)] #, PowerLine(from=1,to=4,K=6)]

powergrid = PowerGrid(nodes,lines)

# TODO: these lines could be avoided: see issue #71
begin
	dimension(myPV)
	symbolsof(myPV)

	dimension(myLoad)
	symbolsof(myLoad)

	#dimension(mySlack)
	#symbolsof(mySlack)
end
systemsize(powergrid)

# using NLsolve: nlsolve, converged
# system_size = systemsize(powergrid)
# ic_guess = ones(system_size)
# rr = RootRhs(rhs(powergrid))
# nl_res = nlsolve(rr, ic_guess)

operationpoint = find_operationpoint(powergrid)

### ILC control algorithm
#= Now putting all things together and defining the ODE
nd = network_dynamics(nodes, lines, powergrid)

begin
	# Define initial values as ic
	# ic = ...
	# Define tspan:
	tspan = (0., num_days * l_day)
	ode_problem = ODEProblem(nd, ic, tspan, compound_pars,
	callback=CallbackSet(PeriodicCallback(HourlyUpdate(), l_hour),
						 PeriodicCallback(DailyUpdate_X, l_day)))
end

@time sol = solve(ode_problem, Rodas4())

for i=1:24*num_days+1
	# get hourly_energy from solver
	hourly_energy[i] = sol((i-1)*3600)[energy_filter[1]]
	# hourly_energy[i,2] = sol1((i-1)*3600)[energy_filter[2]]
	# hourly_energy[i,3] = sol1((i-1)*3600)[energy_filter[3]]
	# hourly_energy[i,4] = sol1((i-1)*3600)[energy_filter[4]]
end

# Apply filtering to get ILC energy for the next days
# ILC_Power is zero for day 1 (day 0, index 1)
ILC_power = zeros(num_days+2,24,N)
# ILC_Power for day 2:
ILC_power[2,:,1] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,1])
# ILC_power[2,:,2] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,2])
# ILC_power[2,:,3] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,3])
# ILC_power[2,:,4] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,4])

for i=2:num_days
	ILC_power[i+1,:,1] = Q*(ILC_power[i,:,1] +  kappa*hourly_energy[(i-1)*24+1:i*24,1])
	# ILC_power[i+1,:,2] = Q*(ILC_power[i,:,2] +  kappa*hourly_energy[(i-1)*24+1:i*24,2])
	# ILC_power[i+1,:,3] = Q*(ILC_power[i,:,3] +  kappa*hourly_energy[(i-1)*24+1:i*24,3])
	# ILC_power[i+1,:,4] = Q*(ILC_power[i,:,4] +  kappa*hourly_energy[(i-1)*24+1:i*24,4])
end
=#

end
