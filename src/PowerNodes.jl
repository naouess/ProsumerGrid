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

    function rhs!(dx, x, e_s, e_d, p, t)
		ϕ = @view x[1]
		ω = @view x[2]
		χ = @view x[3]
		integrated_LI = @view x[4]
		integrated_w = @view x[5]

		u_ILC = p.current_background_power
		F = total_flow(e_s, e_d)
		println("This is the net flow from node PV ", F)

		dϕ = ω[1]
		u_LI = - K_P * ω[1] + χ[1]
		P_gen = (ξ(t) + (u_LI + u_ILC)) * η_gen(t)
		w = 0
		dω = (P_gen - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx[1] = dϕ
		dx[2] = dω
		dx[3] = dχ
		dx[4] = u_LI
		dx[5] = w

		nothing
    end
    ODEVertex(f! = rhs!, dim = 5, mass_matrix = [1, 1, 1, 1, 0], sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
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

	function rhs!(dx, x, e_s, e_d, p, t)
		ϕ = @view x[1]
		ω = @view x[2]
		χ = @view x[3]
		integrated_LI = @view x[4]
		integrated_w = @view x[5]

		u_ILC = p.current_background_power
		F = total_flow(e_s, e_d)
		println("This is the net flow from node load ", F)

		dϕ = ω[1]
		u_LI = - K_P * ω[1] + χ[1]
		w = 0
		P_load = (ξ(t) + (u_LI + u_ILC))/ η_load(t)
		dω = (- P_load - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx[1] = dϕ
		dx[2] = dω
		dx[3] = dχ
		dx[4] = u_LI
		dx[5] = w

		nothing
	end
    ODEVertex(f! = rhs!, dim = 5, mass_matrix = [1, 1, 1, 1, 0], sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
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

    function rhs!(dx, x, e_s, e_d, p, t)
		ϕ = @view x[1]
		ω = @view x[2]
		χ = @view x[3]
		integrated_LI = @view x[4]
		integrated_w = @view x[5]

		u_ILC = p.current_background_power
		F = total_flow(e_s, e_d)
		println("This is the net flow from node slack ", F)

		dϕ = ω[1]
		u_LI = - K_P * ω[1] + χ[1]
		w = 0
		P_gen = (ξ(t) + u_LI + u_ILC) * η_gen(t)
		dχ = (- ω[1] - K_I * χ[1]) * T_inv
		dω = (P_gen - F) * M_inv

		dx[1] = dϕ
		dx[2] = dω
		dx[3] = dχ
		dx[4] = u_LI
		dx[5] = w

		nothing
    end
    ODEVertex(f! = rhs!, dim = 5, mass_matrix = [1, 1, 1, 1, 0], sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
end
