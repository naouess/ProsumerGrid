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
        M
end
PV(; ξ, η_gen, LI, M) = PV(ξ, η_gen, LI, M)
function construct_vertex(par::PV)
    ξ = par.ξ
    η_gen = par.η_gen
    K_P = par.LI.kp
	K_I = par.LI.ki
	T = par.LI.T
    M = par.M

    function rhs!(dx, x, e_s, e_d, p, t)
        ϕ = x[1]
        ω = x[2]
		χ = x[3]
		integrated_LI = x[4]
		integrated_w = x[5]

		u_LI = - K_P * ω + χ
		w = ξ(t) - (u_LI + p.current_background_power)
        P_gen = (u_LI + p.current_background_power) * η_gen(t)
		dχ = (- ω - K_I * χ) * T
		# dϕ = ω
        dω = (P_gen - total_current(e_s, e_d)) * M # F = total_current(e_s, e_d)
        try
            dx[1] = ω
            dx[2] = dω
			dx[3] = dχ
			dx[4] = u_LI
			dx[5] = w
            return nothing
        catch e
            if typeof(e) === UndefVarError
                throw(NodeDynamicsError("you need to provide $(e.var)"))
            else
                throw(e)
            end
        end
    end
    ODEVertex(f! = rhs!, dim = 5, mass_matrix = [1, 1, 1, 1, 1], sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
end

struct Load <: AbstractNode
        ξ
        η_load
        LI
        M
end

Load(; ξ, η_load, LI, M) = Load(ξ, η_load, LI, M)
function construct_vertex(par::Load)
	ξ = par.ξ
    η_load = par.η_load
    K_P = par.LI.kp
	K_I = par.LI.ki
	T = par.LI.T
    M = par.M

	function rhs!(dx, x, e_s, e_d, p, t)
		ϕ = x[1]
		ω = x[2]
		χ = x[3]
		integrated_LI = x[4]
		integrated_w = x[5]

		u_LI = - K_P * ω + χ
		w = ξ(t) - (u_LI + p.current_background_power)
		P_load = (u_LI + p.current_background_power) / η_load(t)
		dχ = (- ω - K_I * χ) * T
		# dϕ = ω
		dω = (- P_load - total_current(e_s, e_d)) * M # F = total_current(e_s, e_d)
		try
			dx[1] = ω
			dx[2] = dω
			dx[3] = dχ
			dx[4] = u_LI
			dx[5] = w
			return nothing
		catch e
			if typeof(e) === UndefVarError
				throw(NodeDynamicsError("you need to provide $(e.var)"))
			else
				throw(e)
			end
		end
	end
	ODEVertex(f! = rhs!, dim = 5, mass_matrix = [1, 1, 1, 1, 1], sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
end

struct Slack <: AbstractNode
        ξ
        η_gen
        LI
        M
end

Slack(; ξ, η_gen, LI, M) = Slack(ξ, η_gen, LI, M)
function construct_vertex(par::Slack)
    ξ = par.ξ
    η_gen = par.η_gen
    K_P = par.LI.kp
	K_I = par.LI.ki
	T = par.LI.T
    M = par.M

    function rhs!(dx, x, e_s, e_d, p, t)
        ϕ = x[1]
        ω = x[2]
		χ = x[3]
		integrated_LI = x[4]
		integrated_w = x[5]

		u_LI = - K_P * ω + χ
		w = ξ(t) - (u_LI + p.current_background_power)
        P_gen = (u_LI + p.current_background_power) * η_gen(t)
		dχ = (- ω - K_I * χ) * T
		# dϕ = ω
        dω = (P_gen - total_current(e_s, e_d)) * M # F = total_current(e_s, e_d)
        try
            dx[1] = ω
            dx[2] = dω
			dx[3] = dχ
			dx[4] = u_LI
			dx[5] = w
            return nothing
        catch e
            if typeof(e) === UndefVarError
                throw(NodeDynamicsError("you need to provide $(e.var)"))
            else
                throw(e)
            end
        end
    end
    ODEVertex(f! = rhs!, dim = 5, mass_matrix = [1, 1, 1, 1, 1], sym = [:ϕ, :ω, :χ, :integrated_LI, :integrated_w])
end
