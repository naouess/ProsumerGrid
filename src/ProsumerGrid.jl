#=
TODO:
[] remove unneeded imports and "using.."-statements
[] correct error in find_operationpoint
=#
module ProsumerGrid

using Parameters
using Pkg
Pkg.instantiate()
using PowerDynamics
import PowerDynamics: dimension, symbolof, construct_vertex

begin
	dir = @__DIR__
	include("$dir/nodes/NodeTypesPD.jl")
	include("$dir/lines/PowerLine.jl")
end

#=
Minimal example to test out the defined nodes and Powerline
No ILC solving integrated yet. u is set as a parameter.
=#
@with_kw mutable struct η
	gen
	load
	"""
	Constructor
	"""
	function η(gen, load)
			new()
	end
end
@with_kw mutable struct u
	ILC
	LI
	"""
	Constructor
	"""
	function u(ILC, LI)
			new()
	end
end
myPV = PV(ξ=1000, η=η(gen=1., load=1.), u=u(ILC=900, LI=900), M=1)
myLoad = Load(ξ=-900, η=η(gen=1., load=1.), u=u(ILC=-900, LI=-900), M=1)
#myBattery = Battery()

nodes = [myPV, myLoad]

lines=[PowerLine(from=1,to=2,K=10)]

powergrid = PowerGrid(nodes,lines)

dimension(myPV)
symbolsof(myPV)

dimension(myLoad)
symbolsof(myLoad)

systemsize(powergrid)

# Here error occurs:
# UndefRefError: access to undefined reference
operationpoint = find_operationpoint(powergrid)

end
