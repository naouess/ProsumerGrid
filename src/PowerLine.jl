abstract type AbstractLine end

begin
    struct PowerLine <: AbstractLine
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
            # If current is flowing away from the source, it is negative at the source.
            source_ϕ = v_s[1]
            destination_ϕ = v_d[1]
            # e = K * sin(destination_ϕ - source_ϕ)
            # Linearized model:
            F = K .* (destination_ϕ .- source_ϕ)
            e .= [K, destination_ϕ, source_ϕ, F]
            nothing
        end
        return StaticEdge(f! = rhs!, dim = 4)
    end
end

### Helper function
function total_flow(x, e_s, e_d)
    # Keeping with the convention of negative sign for outging flow
    net_flow = 0.
    for e in e_s
        net_flow -= e[1] * (x[1] - e[2])
    end
    for e in e_d
        net_flow += e[1] * (e[3] - x[1])
    end
    return net_flow
end
