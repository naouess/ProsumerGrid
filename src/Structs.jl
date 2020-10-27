# TODO:
# Define outer constructors for η for the cases when only gen or load is defined
mutable struct η
	gen
	load
end
η(; gen, load) = η(gen, load)

# TODO: define default values
mutable struct LI
	kp
	ki
	T
end
LI(; kp, ki, T) = LI(kp, ki, T)

mutable struct ILC
	kappa
	mismatch_yesterday
	daily_background_power
	current_background_power
	ilc_nodes
	ilc_cover
	Q
end
ILC(; kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes = 0, ilc_cover = 0, Q) = ILC(kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q)

# export ILC, LI, η
