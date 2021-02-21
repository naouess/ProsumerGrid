abstract type AbstractNode end

## Define node structs
mutable struct PV <: AbstractNode
        ξ
        η_gen
        LI :: LI
		ILC :: ILC
        M
		max_inverter
end
PV(; ξ, η_gen, LI, ILC, M, max_inverter = 0.) = PV(ξ, η_gen, LI, ILC, M, max_inverter)

# TODO finish node definition
mutable struct Wind <: AbstractNode
        ξ
        η_gen
        LI :: LI
		ILC :: ILC
        M
end
Wind(; ξ, η_gen, LI, ILC, M) = Wind(ξ, η_gen, LI, ILC, M)

mutable struct Slack <: AbstractNode
        η_gen
        LI :: LI
		ILC :: ILC
        M
		max
end
Slack(; η_gen, LI, ILC, M, max = 0) = Slack(η_gen, LI, ILC, M, max)

mutable struct Load <: AbstractNode
        ξ
        η_load
        LI :: LI
		ILC :: ILC
        M
end
Load(; ξ, η_load, LI, ILC, M) = Load(ξ, η_load, LI, ILC, M)

mutable struct Battery <: AbstractNode
        η :: η
        LI :: LI
		ILC :: ILC
        M
		C
		max_inverter
		σ
		SOC_min
		SOC_max
end
Battery(; η, LI, ILC, M, C, σ = 0., SOC_min = 0., SOC_max = 1., max_inverter=0.) = Battery(η, LI, ILC, M, C, σ, SOC_min, SOC_max, max_inverter)

mutable struct ThermalStorage <: AbstractNode
        η :: η
        LI :: LI
		ILC :: ILC
        M
		C
		a_loss
		l_ss
		SOC_min
		SOC_max
end
ThermalStorage(; η, LI, ILC, M, C, a_loss = 0., l_ss = 0.,SOC_min = 0., SOC_max = 1.) = ThermalStorage(η, LI, ILC, M, C, a_loss, l_ss)

mutable struct GenericNode <: AbstractNode
		ξ
        η :: η
        LI :: LI
		ILC :: ILC
        M
		C
		max_inverter
		σ
		SOC_min
		SOC_max
end
GenericNode(; ξ, η, LI, ILC, M, C, σ = 0., SOC_min = 0., SOC_max = 1., max_inverter = 0.) = Battery(ξ, η, LI, ILC, M, C, σ, SOC_min, SOC_max, max_inverter)

# export PV, Load, Slack, Battery ...

## Define generic node
function (gn::GenericNode)(dx, x, e_d, p, t)
	u_ILC = gn.ILC.current_background_power
	u_LI = - gn.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = 0
	@inbounds for e in e_d
		F += e[1]
	end

	if x[5] >= gn.SOC_max
		F = min(0., F)
	end
	if x[5] <= gn.SOC_min
		F = max(0., F)
	end

	if  gn.max_inverter != 0. && P_control >= gn.max_inverter
		F = max(- gn.max_inverter, F)
	end

	dx[1] = x[2]
	dx[2] = (gn.ξ(t) + (u_LI + u_ILC) / gn.η_gen(t) + F) / pv.M
	dx[3] = (- x[2] - pv.LI.ki * x[3]) / pv.LI.T
	dx[4] = u_LI
	dx[6] = abs(u_LI)

	if F < 0
		dx[5] = F / gn.η.gen(t) / batt.C / 3600 - gn.σ * x[5]
	else
		dx[5] = F * gn.η.load(t) / batt.C / 3600 - gn.σ * x[5]
	end

	nothing
end

## Define node dynamics depending on type
function (pv::PV)(dx, x, e_d, p, t)
	u_ILC = pv.ILC.current_background_power
	u_LI = - pv.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = 0
	@inbounds for e in e_d
		F += e[1]
	end

	if pv.max_inverter != 0. && P_control >= pv.max_inverter
		F = max(- pv.max_inverter, F)
	end

	dx[1] = x[2]
	dx[2] = (pv.ξ(t) + (u_LI + u_ILC) / pv.η_gen(t) + F) / pv.M
	dx[3] = (- x[2] - pv.LI.ki * x[3]) / pv.LI.T
	dx[4] = u_LI
	dx[5] = abs(u_LI)

	nothing
end

function (load::Load)(dx, x, e_d, p, t)
	u_ILC = load.ILC.current_background_power
	u_LI = - load.LI.kp * x[2] + x[3]

	P_load = u_LI + u_ILC

	F = 0
	@inbounds for e in e_d
		F += e[1]
	end

	dx[1] = x[2]
	dx[2] = (- load.ξ(t) + P_load * load.η_load(t) + F) / load.M
	dx[3] = (- x[2] - load.LI.ki * x[3]) / load.LI.T
	dx[4] = u_LI
	dx[5] = abs(u_LI)
	nothing
end

function (slack::Slack)(dx, x, e_d, p, t)
	u_ILC = slack.ILC.current_background_power
	u_LI = - slack.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = 0
	@inbounds for e in e_d
		F += e[1]
	end

	if slack.max != 0. && P_control >= slack.max
		F = max(-slack.max, F)
	end

	dx[1] = x[2]
	dx[2] = (P_control + F) * slack.M
	dx[3] = (- x[2] - slack.LI.ki * x[3]) / slack.LI.T
	dx[4] = u_LI
	dx[5] = abs(u_LI)

	nothing
end

function (batt::Battery)(dx, x, e_d, p, t)
	u_ILC = batt.ILC.current_background_power
	u_LI = - batt.LI.kp * x[2] + x[3]
	P_control = u_LI + u_ILC

	F = 0
	@inbounds for e in e_d
		F += e[1]
	end

	if x[5] >= batt.SOC_max
		F = min(0., F)
	end
	if x[5] <= batt.SOC_min
		F = max(0., F)
    end

	if batt.max_inverter != 0. && P_control >= batt.max_inverter
		F = max(- batt.max_inverter, F)
	end

	dx[1] = x[2]
	dx[2] = (P_control + F) / batt.M
	dx[3] = (- x[2] - batt.LI.ki * x[3]) / batt.LI.T
	dx[4] = u_LI
	dx[6] = abs(u_LI)

	if F < 0
		dx[5] = - P_control / batt.η.gen(t) / batt.C / 3600 - batt.σ * x[5]
	else
		dx[5] = - P_control * batt.η.load(t) / batt.C / 3600 - batt.σ * x[5]
	end
	nothing
end

function (ts::ThermalStorage)(dx, x, e_d, p, t)
	u_ILC = ts.ILC.current_background_power
	u_LI = - ts.LI.kp * x[2] + x[3]

	P_control = u_LI + u_ILC

	F = 0
	@inbounds for e in e_d
		F += e[1]
	end

	if x[5] <= ts.SOC_min
		F = max(0., F)
	end
	if x[5] >= ts.SOC_max
		F = min(0., F)
	end

	dx[1] = x[2]
	dx[2] = (u_LI + u_ILC + F) / ts.M
	dx[3] = (- x[2] - ts.LI.ki * x[3]) / ts.LI.T
	dx[4] = u_LI
	dx[6] = abs(u_LI)

	if F < 0
		dx[5] =  - P_control / ts.η.gen(t) / ts.C / 3600 - ts.a_loss * (x[5] - ts.l_ss)
	else
		dx[5] =  - P_control * ts.η.load(t) / ts.C / 3600 - ts.a_loss * (x[5] - ts.l_ss)
	end

	nothing
end
