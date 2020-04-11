"""Abstract type for all nodes"""
abstract type AbstractNode end

"""Abstract subtypes. A node can be a bus or a component"""
abstract type AbstractBus <: AbstractNode end
abstract type AbstractComponent <: AbstractNode end

###
using Parameters
# using Base.Meta
# using MacroTools
# ...

struct GenericComponent <: AbstractComponent
    C		# Storage capacity
    η_load	# Load conversion efficiency
    η_gen	# Generation efficiency
    u_load	# Power taken from the grid at node level
    u_gen 	# Power fed into the grid at node level
    du_load	# Ramp
    du_gen 	# Ramp
    l		# State of charge (level)
    dl		# Change in state of charge in time
    ξ		# Provided or demanded energy (depending on sign) => Energy exchange with environment
    w		# Spilled energy or unserved (shedded) load
    v		# Storage losses
	"Inner Constructor for GenericComponent"
	# Defines the default constructor
	function generic_component()
		# evtl. condition to evaluate inputs
		# additional prep work for parameters
		# evtl. any default values on parameters
		# differentiate if some parameters are single values or time series
		new()
	end
end

function main_equation_component(component <: AbstractComponent)
	# Definition of main equation here
	@. component.dl = (# Write rhs of equation here
						) / component.C
	quote

	end
end

function main_equation_bus(bus <: AbstractBus)
	# Definition of main equation here
end

# bounds = [min, max]
function bounds(var, bounds)
	# Bounds on u_gen and u_load
	# Bounds on ramp rates du_gen and du_load
	# Add conditionals
	@assert var in bounds

	# ? How to integrate this ? CallBack-Function ? or using JumP ?
end

### Try with Metaprogramming

parameters = [
	:C,			# Storage capacity
	:η_load,	# Load conversion efficiency
	:η_gen,		# Generation efficiency
	:u_load,	# Power taken from the grid at node level
	:u_gen,		# Power fed into the grid at node level
	:du_load,	# Ramp
	:du_gen, 	# Ramp
	:l,			# State of charge
	:dl,		# Change of SoC
	:ξ,			# Provided or demanded energy (depending on sign) => Energy exchange with environment
	:w,			# Spilled energy or unserved (shedded) load
	:v			# Storage losses
	]

function buildstruct_component(name, parameters)
	struct_def = Expr(
        :struct, false, # why is this second argument needed?
        :($name <: AbstractComponent),
        Expr(:block, parameters...,
            Expr(:(=), # inner constructor
                Expr(:call, name, Expr(:parameters, parameters... )),
                Expr(:call, :new,  parameters...))
    )
end

macro PowerNode(name, parameters)
    # Defnitions
	struct_def_component = buildstruct_component(name, parameters)
	struct_def_bus = buildstruct_bus(name, parameters)
    return
end # macro
