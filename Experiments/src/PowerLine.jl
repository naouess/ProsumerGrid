abstract type AbstractLine end

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
    function edgefunction!(e, v_s, v_d, p, t)
        e[1] = K * sin(v_s[1] - v_d[1])
        nothing
    end
    return StaticEdge(f! = edgefunction!, dim = 1)
end
