# TODO:
# Define outer constructors for η for the cases when only gen or load is defined
mutable struct η
	gen
	load
end
η(; gen, load) = η(gen, load)

mutable struct LI
	kp
	ki
	T_inv
end
LI(; kp, ki, T_inv) = LI(kp, ki, T_inv)

mutable struct ILC
	kappa
	mismatch_yesterday
	daily_background_power
	current_background_power
	ilc_nodes
	ilc_cover
	Q
end
ILC(; kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q) = ILC(kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q)

# export ILC, LI, η

## Not used
# mutable struct compound_pars
# 	ξ
# 	ILC::ILC
# 	LI::LI
# 	M_inv
# 	"""
# 	Constructor
# 	"""
# 	function compound_pars(ξ, ILC, LI, M_inv)
# 			new(ξ, ILC, LI, M_inv)
# 	end
# end
