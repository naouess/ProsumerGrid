

function edgefunction!(e, v_s, v_d, p, t)
    # If current is flowing away from the source, it is negative at the source.
    source_ϕ = v_s[1]
    destination_ϕ = v_d[1]
    # e = K * sin(destination_ϕ - source_ϕ)
    # Linearized model:
    e .= 6. .* sin(source_ϕ - destination_ϕ)
    # e .= [6., destination_ϕ, source_ϕ, flow]
    # println(e)
    nothing
end

line = StaticEdge(f! = edgefunction!, dim = 1)


### Helper function
function total_flow(x, e_s, e_d)
    # Keeping with the convention of negative sign for outging flow
    net_flow = 0.
    for e in e_s
        net_flow += e[1] #* (x - e[2])
    end
    for e in e_d
        net_flow -= e[1] #* (e[3] - x)
    end
    return net_flow
end
