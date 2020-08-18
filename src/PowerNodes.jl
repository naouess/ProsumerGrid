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

		dx1 = @view dx[1]
		dx2 = @view dx[2]
		dx3 = @view dx[3]
		dx4 = @view dx[4]
		dx5 = @view dx[5]

		u_ILC = p.current_background_power
		F = total_flow(ϕ, e_s, e_d)
		# println("This is the net flow from node PV ", F)

		dϕ = ω[1]
		u_LI = - K_P * ω[1] + χ[1]
		P_gen = (ξ(t)[1] - (u_LI + u_ILC)) * η_gen(t)
		# w = ξ(t)[1] - (u_LI + u_ILC)
		w = 0
		dω = (P_gen - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx1 = dϕ
		dx2 = dω
		dx3 = dχ
		dx4 = u_LI
		dx5 = w

		nothing
    end
    ODEVertex(f! = rhs!, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w]) # mass_matrix = [1, 1, 1, 1, 0],
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

		dx1 = @view dx[1]
		dx2 = @view dx[2]
		dx3 = @view dx[3]
		dx4 = @view dx[4]
		dx5 = @view dx[5]

		u_ILC = p.current_background_power
		F = total_flow(ϕ, e_s, e_d)
		# println("This is the net flow from node load ", F)

		dϕ = ω[1]
		u_LI = - K_P * ω[1] + χ[1]
		# w = ξ(t)[1] - (u_LI + u_ILC)
		w = 0
		P_load = (ξ(t)[1] - (u_LI + u_ILC))/ η_load(t)
		dω = (- P_load - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx1 = dϕ
		dx2 = dω
		dx3 = dχ
		dx4 = u_LI
		dx5 = w

		nothing
	end
    ODEVertex(f! = rhs!, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
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

		dx1 = @view dx[1]
		dx2 = @view dx[2]
		dx3 = @view dx[3]
		dx4 = @view dx[4]
		dx5 = @view dx[5]

		u_ILC = p.current_background_power
		F = total_flow(ϕ, e_s, e_d)
		# println("This is the net flow from node Slack ", F)

		dϕ = ω[1]
		u_LI = - K_P * ω[1] + χ[1]
		P_gen = (ξ(t)[1] - (u_LI + u_ILC))* η_gen(t)
		# w = ξ(t)[1] - (u_LI + u_ILC)
		w = 0
		dω = (P_gen - F) * M_inv
		dχ = (- ω[1] - K_I * χ[1]) * T_inv

		dx1 = dϕ
		dx2 = dω
		dx3 = dχ
		dx4 = u_LI
		dx5 = w

		nothing
    end
    ODEVertex(f! = rhs!, dim = 5, sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
end
