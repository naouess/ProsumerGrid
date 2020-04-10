"""Abstract type for all nodes"""
abstract type AbstractNode end

# ? is this additional specification of types needed ?
# ? Maybe reconsider naming of "Bus" and "Component" ?
# Idea behind this: the differentiation is mainly due to the fact that components and buses have different
# equations. So more efficient to seperate both types and have this integrated in
# the type hierarchy.
"""Abstract subtypes. A node can be a bus or a component"""
abstract type AbstractBus <: AbstractNode end
abstract type AbstractComponent <: AbstractNode end

###
# Here import Base packages needed for creating macros
using Parameters
# using Base.Meta
# using MacroTools
# ...

# ? or immutable struct ?
# ? @with_kw in what way does this differ from just accessing the fields by .field_name ?
# ? Use parametric struct to differentiate if some parameters are passed as constants (single values)
# or as time series ?
# Parametric struct => T depends ond whether inputs are constants or time series
mutable struct GenericComponent{T} <: AbstractComponent
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
	# ? dl needs broadcast ?
	@. component.dl = () / component.C

	# use quote end ?
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

function check_input()
	# Here more constraints on input validity check
	#
end

### Try with Metaprogramming (using macro)

# The following is used to declare a struct at parsing time (=> Metaprogramming).
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
        Expr(:block, parameters..., # set all the parmeters as fields in the struct
            Expr(:(=), # define the constructor
                Expr(:call, name, Expr(:parameters, parameters... )),
                Expr(:call, :new,  parameters...))
    )
	# This function returns a piece of code that generate an expression
	# Still unclear how to evaluate this to create the struct at runtime, so that
	# typeof(struct_def) is DataType (CompositeType)

	# ? How to define an inner consructor using quoting => see above! :)
end

# Example of struct definition
name = :Battery
parameters = [:a, :b]
Battery = buildstruct_component(name, parameters)
typeof(Battery) # Type is Expr => still not evaluated

macro PowerNode(name, parameters)
    # Defnitions
	struct_def_component = buildstruct_component(name, parameters)
	struct_def_bus = buildstruct_bus(name, parameters)
    return
end # macro
