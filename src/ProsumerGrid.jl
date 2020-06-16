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

myPV = PV(ξ=800., η_gen=1., LI = u_PV.LI, ILC = u_PV.ILC, M=1.)
myLoad = Load(ξ=-900., η_load = 1., LI = u_load.LI, ILC = u_load.ILC, M=1.)
mySlack = Slack(ξ=100., η_gen=1., LI = u_slack.LI, ILC = u_slack.ILC, M=1.)
slack2 = SlackAlgebraic(U=1.) # To avoid triggering warnings
#Battery = Battery()

nodes = [myPV, myLoad, mySlack, slack2]

lines=[PowerLine(from=1,to=2,K=6), PowerLine(from=2,to=3,K=6), PowerLine(from=1,to=4,K=6)]

powergrid = PowerGrid(nodes,lines)

# TODO: could be avoided.
begin
	dimension(myPV)
	symbolsof(myPV)

	dimension(myLoad)
	symbolsof(myLoad)

	dimension(mySlack)
	symbolsof(mySlack)
end
systemsize(powergrid)

# Error occurs here:
# BoundsError: attempt to access 14-element Array{Float64,1} at index [15]
operationpoint = find_operationpoint(powergrid)
# solution = simulate(LineFault(), powergrid, operationpoint, timespan = (0.0,100.0))
# PowerDynSolve.solve # ?

end
