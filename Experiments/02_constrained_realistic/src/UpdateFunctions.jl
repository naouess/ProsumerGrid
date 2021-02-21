@doc """
    HourlyUpdate(integrator)
PeriodicCallback function acting on the `integrator` that is called every simulation hour (t = 1,2,3...).
"""
function HourlyUpdate(integrator)
	indices = idx_containing(integrator.f, :integrated_LI)

	hour = mod(round(Int, integrator.t/3600.), 24) + 1
	last_hour = mod(hour-2, 24) + 1

	for idx in 1:length(indices)
		integrator.f.f.vertices![idx].f!.ILC.mismatch_yesterday[last_hour] = integrator.u[indices[idx]]
		integrator.u[indices[idx]] = 0.
		integrator.f.f.vertices![idx].f!.ILC.current_background_power = integrator.f.f.vertices![idx].f!.ILC.daily_background_power[hour]
	end

	# indices_w = idx_containing(nd, :curtailed)
	# for idx in 1:length(indices_w)
	# 	integrator.u[indices_w[idx]] = 0.
	# end

	nothing
end

@doc """
    DailyUpdate(integrator)
PeriodicCallback function acting on the `integrator` that is called every simulation day (t = 24,48,..).
"""
function DailyUpdate(integrator)
	num_v = length(integrator.f.f.vertices!)
	for idx in 1:num_v
		integrator.f.f.vertices![idx].f!.ILC.daily_background_power = integrator.f.f.vertices![idx].f!.ILC.Q * (integrator.f.f.vertices![idx].f!.ILC.daily_background_power
		+ integrator.f.f.vertices![idx].f!.ILC.kappa * integrator.f.f.vertices![idx].f!.ILC.mismatch_yesterday)
	end
	nothing
end

function WatchStorage(period)
	deltat = period / 3600.

	function affect!(integrator)
		indices_storage = idx_containing(integrator.f, :level)
		for idx in 1:length(indices_storage)
			integrator.f.f.vertices![4].f!.C = integrator.f.f.vertices![4].f!.C + integrator.u[indices_storage[idx]] * deltat
			integrator.u[indices_storage[idx]] = 0.
		end
	end
	return PeriodicCallback(affect!, period)#, initial_affect= false)
end

# function WatchStorage()
# 	function condition(u, t, integrator)
# 		true
# 	end
# 	function affect!(integrator)
# 		deltat = integrator.t - integrator.tprev
# 		deltat_hour = deltat # / 3600.
#
# 		indices_storage = idx_containing(integrator.f, :level)
# 		for idx in 1:length(indices_storage)
# 			integrator.f.f.vertices![4].f!.C = integrator.f.f.vertices![4].f!.C  + integrator.u[indices_storage[idx]]# * deltat_hour
# 			# integrator.u[indices_storage[idx]] = 0.
# 		end
# 	end
# 	return DiscreteCallback(condition, affect!)
# end
