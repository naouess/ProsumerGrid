using PowerDynamics


# The following import lines are needed to avoid some LoadErrors: there seems to be a
# more proper way to define custom nodes without having to use these imports:
# TODO: Follow up on this (see thread here: https://github.com/JuliaEnergy/PowerDynamics.jl/issues/71)
import Base: @__doc__
import PowerDynamics: AbstractNode, showdefinition, AbstractLine, construct_vertex

# TODO: Check where a broadcast is necessary (parameter and variable types)

@DynamicNode PV(ξ, η_gen, LI, ILC, M) begin
    MassMatrix(; m_u = false, m_int = [true, true])
end begin
    # [all prepratory things that need to be run just once]
    w = ξ - (LI + ILC)
    #@assert ξ * w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    #@assert abs(ξ) - abs(w) >= 0 "ξ always bigger than w"
    #@assert η in [0.0, 1.0]
    # define P as struct?
    P_gen = 0
end [[ϕ, dϕ], [ω, dω]] begin
    # [the actual dynamics equation]
    # [important to set the output variables]
    du = P_gen - (LI + ILC) * η_gen
    F = i
    dϕ = ω
    dω = (P_gen - F) / M
end

@DynamicNode Slack(ξ, η_gen, LI, ILC, M) begin
    MassMatrix(; m_u = false, m_int = [true, true])
end begin
    # [all prepratory things that need to be run just once]
    w = ξ - (LI + ILC)
    #@assert ξ * w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    #@assert abs(ξ) - abs(w) >= 0 "ξ always bigger than w"
    #@assert η in [0.0, 1.0]
    # define P as struct?
    P_gen = 0
end [[ϕ, dϕ], [ω, dω]] begin
    # [the actual dynamics equation]
    # [important to set the output variables]
    du = P_gen - (LI + ILC) * η_gen
    F = i
    dϕ = ω
    dω = (P_gen - F) / M
end

@DynamicNode Load(ξ, η_load, LI, ILC, M) begin
    MassMatrix(; m_u = false, m_int = [true, true])
end begin
    w = ξ - (LI + ILC)
    #@assert ξ * w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    #@assert abs(ξ) - abs(w) >= 0 "ξ always bigger than w"
    #@assert η in [0.0, 1.0]
    # define P as struct?
    P_load = 0
end [[ϕ, dϕ], [ω, dω]] begin
    du = P_load - (LI + ILC) / η_load
    F = i # since i is the variable for the flows as defined in NodeMacro.jl
    dϕ = ω
    dω = (- P_load - F) / M
end

# TODO: Recheck definition of battery node
@DynamicNode Battery(η, M, C) begin
    MassMatrix(; m_u = false, m_int = [true, true])
end begin
    @assert η in [0.0, 1.0]
    P.load = []
    P.gen = []
end [[ϕ, dϕ], [ω, dω]] begin
    # ? Differentiate between charging and discharging states of battery through a
    # conditional on the sign of F? (both states don't occur simultaneously)
    # P.gen = (u.LI + u.ILC) * η.gen
    # P.load = (u.LI + u.ILC) / η.load
    # ? Can du be defined like this?
    du = 0
    dl = (η.load * P.load - P.gen / η.gen) / C
    F = i # since i is the variable for the flows as defined in NodeMacro.jl
    dϕ = ω
    dω = (P.gen - P.load - F) / M
end
