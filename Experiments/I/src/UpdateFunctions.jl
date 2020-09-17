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
	# The current hour is extracted using the following line:
	hour = mod(round(Int, integrator.t/3600.), 24) + 1
	# getting the last hour of the day
	last_hour = mod(hour-2, 24) + 1

	# the following lines: last hour letzte Stunde auf Integral der Stunde davor setzen
	for j = 1:N
		integrator.p[1].mismatch_yesterday[last_hour, j] = integrator.u[4*j]
		integrator.u[4*j] = 0.
		integrator.p[1].current_background_power[j] = integrator.p[1].daily_background_power[hour, j]
	end
	# Why is the abs. value needed ?
	# integrator.u[power_abs_idx] .= 0.
	nothing
end

function DailyUpdate(integrator)
	integrator.p[1].daily_background_power = integrator.p[1].Q * (integrator.p[1].daily_background_power
	+ integrator.p[1].kappa * integrator.p[1].mismatch_yesterday)
	nothing
end
