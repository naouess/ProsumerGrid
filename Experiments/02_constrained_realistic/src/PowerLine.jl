abstract type AbstractLine end

mutable struct PowerLine <: AbstractLine
    from
    to
    K
    PowerLine(; from, to, K) = new(from, to, K)
end

function (edge::PowerLine)(e, v_s, v_d, p, t)
    e[1] = edge.K * sin(v_s[1] - v_d[1])
    nothing
end

# export PowerLine
