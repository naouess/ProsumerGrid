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

	pars = view(integrator.p[1], 1:4)
	j = 1

	# the following lines: last hour letzte Stunde auf Integral der Stunde davor setzen
	for param in pars
		# println("The mismatch for ", last_hour, " at node ", j, " is: ", integrator.p[1][j].ILC.mismatch_yesterday[last_hour])
		# println("The integrator.u[4j] at ", last_hour, " at node ", j, " is: ", integrator.u[4*j])
		# println(integrator.p)
		println("I'm at node ", j)
		println(param.ILC.mismatch_yesterday[last_hour], " and ", integrator.u[4*j])
		param.ILC.mismatch_yesterday[last_hour] = integrator.u[4*j] #???????
		println(param.ILC.mismatch_yesterday[last_hour])
		integrator.u[4*j] = 0.
		param.ILC.current_background_power = param.ILC.daily_background_power[hour]
		j += 1
	end

	# println("The mismatch for ", last_hour, " at node 1 is: ", integrator.p[1][1].ILC.mismatch_yesterday[last_hour])
	# println("The integrator.u[4] at ", last_hour, " at node 1 is: ", integrator.u[4])
	# integrator.p[1][1].ILC.mismatch_yesterday[last_hour] = integrator.u[4]
	# integrator.u[4] = 0.
	# integrator.p[1][1].ILC.current_background_power = integrator.p[1][1].ILC.daily_background_power[hour]
	#
	# println("The mismatch for ", last_hour, " at node 2 is: ", integrator.p[1][2].ILC.mismatch_yesterday[last_hour])
	# println("The integrator.u[8] at ", last_hour, " at node 2 is: ", integrator.u[8])
	# integrator.p[1][2].ILC.mismatch_yesterday[last_hour] = integrator.u[8]
	# integrator.u[8] = 0.
	# integrator.p[1][2].ILC.current_background_power = integrator.p[1][2].ILC.daily_background_power[hour]
	#
	# println("The mismatch for ", last_hour, " at node 3 is: ", integrator.p[1][3].ILC.mismatch_yesterday[last_hour])
	# println("The integrator.u[12] at ", last_hour, " at node 3 is: ", integrator.u[12])
	# integrator.p[1][3].ILC.mismatch_yesterday[last_hour] = integrator.u[12]
	# integrator.u[12] = 0.
	# integrator.p[1][3].ILC.current_background_power = integrator.p[1][3].ILC.daily_background_power[hour]
	#
	# println("The mismatch for ", last_hour, " at node 4 is: ", integrator.p[1][4].ILC.mismatch_yesterday[last_hour])
	# println("The integrator.u[16] at ", last_hour, " at node 4 is: ", integrator.u[16])
	# integrator.p[1][4].ILC.mismatch_yesterday[last_hour] = integrator.u[16]
	# integrator.u[16] = 0.
	# integrator.p[1][4].ILC.current_background_power = integrator.p[1][4].ILC.daily_background_power[hour]

	nothing
end

function DailyUpdate(integrator)
	for j = 1:N
		integrator.p[1][j].ILC.daily_background_power = integrator.p[1][j].ILC.Q * (integrator.p[1][j].ILC.daily_background_power
		+ integrator.p[1][j].ILC.kappa * integrator.p[1][j].ILC.mismatch_yesterday)
	end
	println("I did the daily update")
	nothing
end
