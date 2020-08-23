using NetworkDynamics, Parameters
begin
	dir = @__DIR__
	include("$dir/PowerLine.jl")
	include("$dir/structs.jl")
end
abstract type AbstractNode end

struct PV <: AbstractNode
        ξ
        η_gen
        LI
        M_inv
end
PV(; ξ, η_gen, LI, M_inv) = PV(ξ, η_gen, LI, M_inv)
function construct_vertex(par::PV)
    ξ = par.ξ
    η_gen = par.η_gen
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
    M_inv = par.M_inv

    function PVfunction!(dx, x, e_s, e_d, p, t)
		ϕ = view(x, 1)
		ω = view(x, 2)
		χ = view(x, 3)

		u_ILC = p.current_background_power
		F = total_flow(ϕ, e_s, e_d)

		u_LI = - K_P * ω[1] + χ[1]
		w = ξ(t)[1] - (u_LI + u_ILC)
		P_gen = w * η_gen(t)
		dω = (P_gen - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx[1] = ω[1]
		dx[2] = dω
		dx[3] = dχ
		dx[4] = u_LI

		nothing
    end
    ODEVertex(f! = rhs!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end

struct Load <: AbstractNode
        ξ
        η_load
        LI
        M_inv
end

Load(; ξ, η_load, LI, M_inv) = Load(ξ, η_load, LI, M_inv)
function construct_vertex(par::Load)
	ξ = par.ξ
    η_load = par.η_load
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
	M_inv = par.M_inv

	function Loadfunction!(dx, x, e_s, e_d, p, t)
		ϕ = view(x, 1)
		ω = view(x, 2)
		χ = view(x, 3)

		u_ILC = p.current_background_power
		F = total_flow(ϕ, e_s, e_d)

		u_LI = - K_P * ω[1] + χ[1]
		w = ξ(t)[1] - (u_LI + u_ILC)
		P_load = w / η_load(t)
		dω = (-P_load- F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx[1] = ω[1]
		dx[2] = dω
		dx[3] = dχ
		dx[4] = u_LI

		nothing
    end
    ODEVertex(f! = rhs!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end

struct Slack <: AbstractNode
        ξ
        η_gen
        LI
        M_inv
end

Slack(; ξ, η_gen, LI, M_inv) = Slack(ξ, η_gen, LI, M_inv)
function construct_vertex(par::Slack)
    ξ = par.ξ
    η_gen = par.η_gen
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
	M_inv = par.M_inv

    function Slackfunction!(dx, x, e_s, e_d, p, t)
		ϕ = view(x, 1)
		ω = view(x, 2)
		χ = view(x, 3)

		u_ILC = p.current_background_power
		F = total_flow(ϕ, e_s, e_d)

		u_LI = - K_P * ω[1] + χ[1]
		w = ξ(t)[1] - (u_LI + u_ILC)
		P_gen = w * η_gen(t)
		dω = (P_gen - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx[1] = ω[1]
		dx[2] = dω
		dx[3] = dχ
		dx[4] = u_LI

		nothing
    end
    ODEVertex(f! = rhs!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end
