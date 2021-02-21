# module ProsumerGrid

## TODO
# (    )  Add dependencies: using ProsumerGrid
# (done) generalize p.current_background_power[]:
# (done) add w component
# (done) try with mass_matrix
# (done) add battery model
# (done) add efficiencies
# (done) add bounds
# (done) adjust i0
# (done) adjust how hourly energy is calculated
# (    ) add bus node
# (    ) make load node uncontrollable
# (done) try with all LI_parameters the same
# (    ) investigate effect of M
# (    ) add generic node model for experimentations :)

using NetworkDynamics
using LightGraphs
using DifferentialEquations

include("Utilities.jl")
include("Structs.jl")
include("PowerLine.jl")
include("PowerNodes.jl")
include("UpdateFunctions.jl")

# export

# end # module
