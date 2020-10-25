abstract type AbstractNode end

## Define node structs
mutable struct PV <: AbstractNode
        ξ
        η_gen
        LI :: LI
		ILC :: ILC
        M_inv
end
PV(; ξ, η_gen, LI, ILC, M_inv) = PV(ξ, η_gen, LI, ILC, M_inv)

mutable struct Wind <: AbstractNode
        ξ
        η_gen
        LI :: LI
		ILC :: ILC
        M_inv
end
Wind(; ξ, η_gen, LI, ILC, M_inv) = Wind(ξ, η_gen, LI, ILC, M_inv)

mutable struct Slack <: AbstractNode
        η_gen
        LI :: LI
		ILC :: ILC
        M_inv
end
Slack(; η_gen, LI, ILC, M_inv) = Slack(η_gen, LI, ILC, M_inv)

mutable struct Load <: AbstractNode
        ξ
        η_load
        LI :: LI
		ILC :: ILC
        M_inv
end
Load(; ξ, η_load, LI, ILC, M_inv) = Load(ξ, η_load, LI, ILC, M_inv)

mutable struct Battery <: AbstractNode
        η :: η
        LI :: LI
		ILC :: ILC
        M_inv
		C
		σ
end
Battery(; η, LI, ILC, M_inv, C, σ = 0) = Battery(η, LI, ILC, M_inv, C, σ)

mutable struct ThermalStorage <: AbstractNode
        η :: η
        LI :: LI
		ILC :: ILC
        M_inv
		C
		a_loss
		l_ss
end
ThermalStorage(; η, LI, ILC, M_inv, C, a_loss, l_ss) = ThermalStorage(η, LI, ILC, M_inv, C, a_loss, l_ss)

# export PV, Load, Slack, Battery ...

## Define node dynamics depending on type
function (pv::PV)(dx, x, e_s, e_d, p, t)
	u_ILC = pv.ILC.current_background_power
	u_LI = - pv.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = total_flow(e_s, e_d)

	if P_control >= 65. / 50.
		F = max(-65/50, F)
	end

	dx[1] = x[2]
	dx[2] = (pv.ξ(t) + (u_LI + u_ILC) / pv.η_gen(t) + F) * pv.M_inv
	dx[3] = (- x[2] - pv.LI.ki * x[3]) * pv.LI.T_inv
	dx[4] = u_LI
	dx[5] = abs(u_LI)

	nothing
end

function (load::Load)(dx, x, e_s, e_d, p, t)
	u_ILC = load.ILC.current_background_power
	u_LI = - load.LI.kp * x[2] + x[3]

	P_load = u_LI + u_ILC

	F = total_flow(e_s, e_d)

	dx[1] = x[2]
	dx[2] = (- load.ξ(t) + P_load * load.η_load(t) + F) * load.M_inv
	dx[3] = (- x[2] - load.LI.ki * x[3]) * load.LI.T_inv
	dx[4] = u_LI
	dx[5] = abs(u_LI)
	nothing
end

function (slack::Slack)(dx, x, e_s, e_d, p, t)
	u_ILC = slack.ILC.current_background_power
	u_LI = - slack.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = total_flow(e_s, e_d)

	if P_control >= 27/ 50
		F = max(-27/50, F)
	end

	dx[1] = x[2]
	dx[2] = (P_control + F) * slack.M_inv
	dx[3] = (- x[2] - slack.LI.ki * x[3]) * slack.LI.T_inv
	dx[4] = u_LI
	dx[5] = abs(u_LI)

	nothing
end

function (batt::Battery)(dx, x, e_s, e_d, p, t)
	u_ILC = batt.ILC.current_background_power
	u_LI = - batt.LI.kp * x[2] + x[3]
	P_control = u_LI + u_ILC

	F = total_flow(e_s, e_d)

	if x[5] >= 0.9
		F = min(0., F)
	end
	if x[5] <= 0.4
		F = max(0., F)
    end

	if P_control >= 49. ./ 50.
		F = max(- 49. ./ 50., F)
	end

	dx[1] = x[2]
	dx[2] = (P_control + F) * batt.M_inv
	dx[3] = (- x[2] - batt.LI.ki * x[3]) * batt.LI.T_inv
	dx[4] = u_LI
	dx[6] = abs(u_LI)

	if F < 0
		dx[5] = - P_control / batt.η.gen(t) / 3600 - batt.σ * x[5]
	else
		dx[5] = - P_control * batt.η.load(t) / 3600 - batt.σ * x[5]
	end
	nothing
end

function (ts::ThermalStorage)(dx, x, e_s, e_d, p, t)
	u_ILC = ts.ILC.current_background_power
	u_LI = - ts.LI.kp * x[2] + x[3]

	F = total_flow(e_s, e_d)
	if x[5] <= 0.
		F = max(0., F)
	end
	if x[5] >= 1.
		F = min(0., F)
	end

	dx[1] = x[2]
	dx[2] = (u_LI + u_ILC + F) * ts.M_inv
	dx[3] = (- x[2] - ts.LI.ki * x[3]) * ts.LI.T_inv
	dx[4] = u_LI
	dx[6] = abs(u_LI)
	if F < 0
		dx[5] =  -(u_LI + u_ILC) / ts.η.gen(t) / ts.C / 3600 # - a_loss * (x[5] - l_ss) # 900
	else
		dx[5] =  -(u_LI + u_ILC) * ts.η.load(t) / ts.C / 3600 # - a_loss * (x[5] - l_ss) # 900
	end

	nothing
end

function (wind::Wind)(dx, x, e_s, e_d, p, t)
	u_ILC = wind.ILC.current_background_power
	u_LI = - wind.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = total_flow(e_s, e_d)
	w = wind.ξ(t) + (u_LI + u_ILC) / wind.η_gen(t)

	dx[1] = x[2]
	dx[2] = (u_LI + u_ILC + F) * wind.M_inv
	dx[3] = (- x[2] - wind.LI.ki * x[3]) * wind.LI.T_inv
	dx[4] = u_LI
	dx[5] = w

	nothing
end

## Define constructors for each node type
function constructor(f::PV)
	# @assert
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_abs_LI])
end
function constructor(f::Wind)
	# @assert
	return ODEVertex(f! = f, dim = 6, sym = [:ϕ, :ω, :χ, :integrated_LI, :curtailed, :integrated_abs_LI])
end
function constructor(f::Load)
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_abs_LI])
end
function constructor(f::Slack)
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_abs_LI])
end

function constructor(f::Battery)
	return ODEVertex(f! = f, dim = 6, sym = [:ϕ, :ω, :χ, :integrated_LI, :level, :integrated_abs_LI])
end

function constructor(f::ThermalStorage)
	return ODEVertex(f! = f, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :level, :integrated_abs_LI])
end
