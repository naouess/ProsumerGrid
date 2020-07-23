using PowerDynamics
using NetworkDynamics

#=
# The following code defines the PowerLine with the @PowerLine macro
begin
	dir = @__DIR__
	include("$dir/PowerLineMacro.jl")
end
import PowerDynamics: AbstractLine

@PowerLine PowerLine(from, to, K) begin
    F_ij = K * sin(destination_ϕ - source_ϕ)
    F_ij_vector = [F_ij, F_ij]
end
=#
begin
    @__doc__ struct PowerLine <: AbstractLine
            from
            to
            K
            PowerLine(; from, to, K) = new(from, to, K)
        end
    function construct_edge(par::PowerLine)
        from = par.from
        to = par.to
        K = par.K
        function rhs!(e, v_s, v_d, p, t)
            source_ϕ = v_s[3]
            destination_ϕ = v_d[3]
            F_ij = K * sin(destination_ϕ - source_ϕ)
            F_ij_vector = [F_ij, F_ij]
            e[1] = F_ij_vector[1]
            e[2] = 0
            e[3] = F_ij_vector[2]
            e[4] = 0
        end
        return StaticEdge(f! = rhs!, dim = 4)
    end
end

export PowerLine
