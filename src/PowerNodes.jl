begin
	dir = @__DIR__
	include("$dir/Structs.jl")
end

## Helper function
@inline function total_flow!(net_flow, e_s, e_d, x)
    @inbounds for e in e_s
        net_flow -= e[1] * x
    end
    @inbounds for e in e_d
        net_flow += e[1] * x
    end
	nothing
end

## Defining construct_vertex function to describe node dynamics
function construct_vertex(par::PV)
    ξ = par.ξ
    η_gen = par.η_gen
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
    M_inv = par.M_inv

    function PVfunction!(dx, x, e_s, e_d, p, t)
		u_ILC = p.current_background_power[1]

		u_LI = - K_P * x[2] + x[3]

		dx[1] = x[2]
		dx[2] = (- ξ(t) + u_LI + u_ILC) * M_inv
		total_flow!(dx[2], e_s, e_d, M_inv)
		dx[3] = (- x[2] - K_I * x[3]) * T_inv
		dx[4] = u_LI

		nothing
	end
    ODEVertex(f! = PVfunction!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end

function construct_vertex(par::Load)
	ξ = par.ξ
    η_load = par.η_load
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
	M_inv = par.M_inv

	function Loadfunction!(dx, x, e_s, e_d, p, t)
		u_ILC = p.current_background_power[2]

		u_LI = - K_P * x[2] + x[3]

		dx[1] = x[2]
		dx[2] = (- ξ(t) + u_LI + u_ILC) * M_inv
		total_flow!(dx[2], e_s, e_d, M_inv)
		dx[3] = (- x[2] - K_I * x[3]) * T_inv
		dx[4] = u_LI

		nothing
    end
    ODEVertex(f! = Loadfunction!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end

function construct_vertex(par::Load2)
	ξ = par.ξ
    η_load = par.η_load
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
	M_inv = par.M_inv

	function Loadfunction2!(dx, x, e_s, e_d, p, t)
		u_ILC = p.current_background_power[4]

		u_LI = - K_P * x[2] + x[3]

		dx[1] = x[2]
		dx[2] = (- ξ(t) + u_LI + u_ILC) * M_inv
		total_flow!(dx[2], e_s, e_d, M_inv)
		dx[3] = (- x[2] - K_I * x[3]) * T_inv
		dx[4] = u_LI

		nothing
    end
    ODEVertex(f! = Loadfunction2!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end

function construct_vertex(par::Slack)
    ξ = par.ξ
    η_gen = par.η_gen
    K_P = par.LI.kp
	K_I = par.LI.ki
	T_inv = par.LI.T_inv
	M_inv = par.M_inv

    function Slackfunction!(dx, x, e_s, e_d, p, t)
		u_ILC = p.current_background_power[3]

		u_LI = - K_P * x[2] + x[3]

		dx[1] = x[2]
		dx[2] = (- ξ(t) + u_LI + u_ILC) * M_inv
		total_flow!(dx[2], e_s, e_d, M_inv)
		dx[3] = (- x[2] - K_I * x[3]) * T_inv
		dx[4] = u_LI

		nothing
	end
    ODEVertex(f! = Slackfunction!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
end
