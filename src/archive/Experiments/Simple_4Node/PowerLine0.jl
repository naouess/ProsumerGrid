function edgefunction!(e, v_s, v_d, p, t)
    e[1] = 6. * sin(v_s[1] - v_d[1])
    nothing
end

line = StaticEdge(f! = edgefunction!, dim = 1)
