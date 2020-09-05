
@with_kw mutable struct η
	gen
	load
	"""
	Constructor
	"""
	function η(gen, load)
			new(gen, load)
	end
end
# TODO: Define outer constructors for η for the cases when only gen or load is defined

@with_kw mutable struct LI
	kp
	ki
	T_inv
	"""
	Constructor
	"""
	function LI(kp, ki, T_inv)
			new(kp, ki, T_inv)
	end
end

@with_kw mutable struct ILC
	kappa
	mismatch_yesterday
	daily_background_power
	current_background_power
	ilc_nodes
	ilc_cover
	Q
	"""
	Constructor
	"""
	function ILC(kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q)
			new(kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q)
	end
end

@with_kw mutable struct compound_pars
	ξ
	ILC::ILC
	LI::LI
	M_inv
	"""
	Constructor
	"""
	function compound_pars(ξ, ILC, LI, M_inv)
			new(ξ, ILC, LI, M_inv)
	end
end


## Define node types as structs:
abstract type AbstractNode end

struct PV <: AbstractNode
        ξ
        η_gen
        LI
        M_inv
end
PV(; ξ, η_gen, LI, M_inv) = PV(ξ, η_gen, LI, M_inv)

struct Slack <: AbstractNode
        ξ
        η_gen
        LI
        M_inv
end

Slack(; ξ, η_gen, LI, M_inv) = Slack(ξ, η_gen, LI, M_inv)

struct Load <: AbstractNode
        ξ
        η_load
        LI
        M_inv
end

Load(; ξ, η_load, LI, M_inv) = Load(ξ, η_load, LI, M_inv)

struct Load2 <: AbstractNode
        ξ
        η_load
        LI
        M_inv
end

Load2(; ξ, η_load, LI, M_inv) = Load2(ξ, η_load, LI, M_inv)
