using PowerDynamics

begin
	dir = @__DIR__
	include("$dir/PowerLineMacro.jl")
end
import PowerDynamics: AbstractLine

@PowerLine PowerLine(from, to, K) begin
    F_ij = K * sin(destination_ω - source_ω)
    F_ij_vector = [F_ij, F_ij]
end

export PowerLine
