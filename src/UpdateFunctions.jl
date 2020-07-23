# Update functions
@doc """
    HourlyUpdate()
Store the integrated control power in memory.
"""
struct HourlyUpdate
	integrated_control_power_history
	HourlyUpdate() = new([])
end

@doc """
    HourlyUpdate(integrator)
PeriodicCallback function acting on the `integrator` that is called every simulation hour (t = 1,2,3...).
"""
function (hu::HourlyUpdate)(integrator)
	# Here the current hour is extracted using the following line:
	hour = mod(round(Int, integrator.t/3600.), 24) + 1
	# Needed for indexing (Eine Stunde geht von 0 bis 1, also letzteres nehmen)
	last_hour = mod(hour-2, 24) + 1

	# the following lines: last hour letzte Stunde auf Integral der Stunde davor setzen
	for j = 1:N
		integrator.p[1][j].mismatch_yesterday[last_hour] = integrator.u[4*j]
		integrator.u[4*j] = 0.
		integrator.p[1][j].current_background_power = integrator.p[1][j].daily_background_power[hour]
	end
	# Why is the abs. value needed ?
	# integrator.u[power_abs_idx] .= 0.
	nothing
end

function DailyUpdate_X(integrator)
	for j = 1:N
		integrator.p[1][j].daily_background_power .= integrator.p[1][j].Q * (integrator.p[1][j].daily_background_power
		.+ integrator.p[1][j].kappa * integrator.p[1][j].mismatch_yesterday)
	end
	nothing
end
