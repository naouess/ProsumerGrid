#=
TODO:
[̌] Check where a broadcast is necessary
[] Define MassMatrix correctly: u has no mass, whereas ω does.
=#

using PowerDynamics

# The following import lines are needed to avoid some LoadErrors: there seems to be a
# more proper way to define custom nodes without having to use these imports:
# Follow up on this (see thread here: https://github.com/JuliaEnergy/PowerDynamics.jl/issues/71)
import Base: @__doc__
import PowerDynamics: AbstractNode
import PowerDynamics: showdefinition

@DynamicNode PV(ξ, η, u, M) begin
    #MassMatrix() ???
end begin
    # [all prepratory things that need to be run just once]
    @. w = ξ - (u.LI + u.ILC)
    @assert ξ .* w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    @assert abs(ξ) .- abs(w) >= 0 "ξ always bigger than w"
    @assert η in [0.0, 1.0]
    # define P as struct?
    P.gen = []
end [[ω, dω]] begin
    # [the actual dynamics equation]
    # [important to set the output variables]
    du = P.gen - (u.LI + u.ILC) * η.gen
    F = i
    dω = (P.gen - F) / M
end

@DynamicNode Load(ξ, η, u, M) begin
    #MassMatrix()???
end begin
    # [all prepratory things that need to be run just once]
    @. w = ξ - (u.LI + u.ILC)
    @assert ξ .* w >= 0 "ξ and w must have the same sign, according to the defined sign convention"
    @assert abs(ξ) .- abs(w) >= 0 "ξ always bigger than w"
    @assert η in [0.0, 1.0]
    # define P as struct?
    P.load = []
end [[ω, dω]] begin
    # [the actual dynamics equation]
    # [important to set the output variables]
    du = P.load - (u.LI + u.ILC) / η.load
    F = i # since i is the variable for the flows as defined in NodeMacro.jl
    dω = (P.gen - F) / M
end

@DynamicNode Battery(ξ, η, u, M) begin
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

end
