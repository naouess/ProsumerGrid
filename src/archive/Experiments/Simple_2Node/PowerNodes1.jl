using NetworkDynamics
begin
	dir = @__DIR__
	include("$dir/PowerLine1.jl")
	include("$dir/Structs.jl")
end

function rhs!(dx, x, e_s, e_d, p, t)
	ϕ = view(x, 1)
	ω = view(x, 2)
	χ = view(x, 3)

	u_ILC = p.ILC.current_background_power
	F = - total_flow(e_s, e_d)

	u_LI = - p.LI.kp * ω[1] + χ[1]

	dx[1] = ω[1]
	dx[2] = (- p.ξ(t) + (u_LI + u_ILC) + F) * p.M_inv
	dx[3] = (- ω[1] - p.LI.ki * χ[1]) * p.LI.T_inv
	dx[4] = u_LI

	nothing
end

node = ODEVertex(f! = rhs!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
