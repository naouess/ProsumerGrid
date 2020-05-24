# PD
using PowerDynamics

# The following lines are needed to avoid some LoadErrors: there seems to be a more proper
# way to define custom nodes without having to use these imports (see thread here:
# https://github.com/JuliaEnergy/PowerDynamics.jl/issues/71)

import Base: @__doc__
import PowerDynamics: AbstractNode
import PowerDynamics: showdefinition

#=
TODO:
[̌] Check where a broadcast is necessary
[] Define MassMatrix correctly: u has no mass, whereas ω does
=#
@DynamicNode PV(ξ, η, u, M) begin
    #MassMatrix() ???
end begin
    # [all prepratory things that need to be run just once]
    @. w = ξ - (u.LI + u.ILC)
    @assert ξ .* w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    @assert abs(ξ) .- abs(w) >= 0 "ξ always bigger than w"
    @assert η in [0.0, 1.0]
    P.gen = []
end [[ω, dω]] begin
    # [the actual dynamics equation]
    # [important to set the output variables]
    du = P.gen - (u.LI + u.ILC) * η.gen
    F = i
    dω = (P.gen - F) / M
end
# export PV

@DynamicNode Load(ξ, η, u, M) begin
    #MassMatrix()???
end begin
    # [all prepratory things that need to be run just once]
    @. w = ξ - (u.LI + u.ILC)
    @assert ξ .* w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    @assert abs(ξ) .- abs(w) >= 0 "ξ always bigger than w"
    @assert η in [0.0, 1.0]
    P.load = []
end [[ω, dω]] begin
    # [the actual dynamics equation]
    # [important to set the output variables]
    du = P.load - (u.LI + u.ILC) * η.gen
    F = i
    dω = (P.gen - F) / M
end
