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
            source_ϕ = v_s[1]
            destination_ϕ = v_d[1]
            # e = K * sin(destination_ϕ - source_ϕ)
            # Linear model:
            e = K * (destination_ϕ - source_ϕ)
        end
        return StaticEdge(f! = rhs!, dim = 1)
    end
end

### Helper function
function total_current(e_s, e_d)
    # Keeping with the convention of negative sign for outging current
    net_flow = 0.
    for e in e_s
        net_flow -= e[1]
    end
    for e in e_d
        net_flow += e[1]
    end
    net_flow
end
