using MacroTools

import Base: @__doc__

macro PowerLine(typedef, prep, func_body)
    return create(typedef, prep, func_body)
end

macro PowerLine(typedef, func_body)
    return create(typedef, nothing, func_body)
end

function create_line(name, prep, parameters, func_body)
    struct_exp = create_struct(name, parameters)

    ce_call = :(construct_edge(par::$(name)))
    extracted_parameters = map(sym -> :( $sym = par.$sym ), parameters)
    cl_body = quote end
    append!(cl_body.args, extracted_parameters)
    if (prep !== nothing)
        append!(cl_body.args, prep.args)
    end
    cl_function = Expr(:function, ce_call, cl_body)

    rhscall = :(rhs!(e,v_s,v_d,p,t))
    rhsbody = quote end
    rhsbody.args[1] = func_body.args[1]
    append!(rhsbody.args, [:(source_ω = v_s[1])])
    append!(rhsbody.args, [:(destination_ω = v_d[1])])
    append!(rhsbody.args, func_body.args)

    es = [:(e[1] = F_ij_vector[1])]
    ed = [:(e[2] = F_ij_vector[2])]

    append!(rhsbody.args, [es; ed])

    rhs_function_exp = Expr(:function, rhscall, rhsbody)
    edge_exp = :(return StaticEdge(f! = rhs!, dim = 2))
    append!(cl_function.args[2].args, [rhs_function_exp, edge_exp])

    ret = quote
        @__doc__ $(struct_exp)
        $(cl_function)
    end
    return ret
end

function create_struct(name, parameters)
    struct_def = Expr(
        :struct, false,
        :($name <: AbstractLine),
        Expr(:block, parameters..., # set all the parmeters as fields in the struct
            Expr(:(=), # define the constructor
                Expr(:call, name, Expr(:parameters, parameters... )),
                Expr(:call, :new,  parameters...)
            )
        )
    )
end

function create_showdefinition(exp, name)
    mainexstr = "$(copy(exp)|>rmlines|> MacroTools.striplines)"
    return :(showdefinition(io::IO, ::Type{$name}) = println(io, $mainexstr))
end

function create(typedef, prep, func_body)
    @capture(typedef, name_(parameters__))
    mainex = create_line(name, prep, parameters, func_body)
    showex = create_showdefinition(mainex, name)
    append!(mainex.args, [showex])
    return esc(mainex)
end
