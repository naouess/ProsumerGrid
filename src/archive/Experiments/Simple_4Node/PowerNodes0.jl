begin
	dir = @__DIR__
	include("$dir/Structs.jl")
end

### Helper function
@inline function total_flow!(net_flow, e_s, e_d, p)
    @inbounds for e in e_s
        net_flow -= e[1] * p.M_inv
    end
    @inbounds for e in e_d
        net_flow += e[1] * p.M_inv
    end
	# net_flow = net_flow * p.M_inv
	nothing
end

function vertexfunction!(dx, x, e_s, e_d, p, t)
	u_ILC = p.ILC.current_background_power

	dx[1] = x[2]
	dx[2] = (- p.ξ(t) + dx[4] + u_ILC) * p.M_inv
	total_flow!(dx[2], e_s, e_d, p)
	dx[3] = (- x[2] - p.LI.ki * x[3]) * p.LI.T_inv
	dx[4] = - p.LI.kp * x[2] + x[3]

	nothing
end

node = ODEVertex(f! = vertexfunction!, dim = 4, sym = [:ϕ, :ω, :χ, :integrated_LI])
