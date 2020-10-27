function HourlyUpdate(integrator)
	indices = idx_containing(integrator.f, :integrated_LI)

	hour = mod(round(Int, integrator.t/3600.), 24) + 1
	last_hour = mod(hour-2, 24) + 1

	for idx in 1:length(indices)
		integrator.f.f.vertices![idx].f!.ILC.mismatch_yesterday[last_hour] = integrator.u[indices[idx]]
		integrator.u[indices[idx]] = 0.
		integrator.f.f.vertices![idx].f!.ILC.current_background_power = integrator.f.f.vertices![idx].f!.ILC.daily_background_power[hour]
	end

	nothing
end

function DailyUpdate(integrator)
	num_v = length(integrator.f.f.vertices!)
	for idx in 1:num_v
		integrator.f.f.vertices![idx].f!.ILC.daily_background_power = integrator.f.f.vertices![idx].f!.ILC.Q * (integrator.f.f.vertices![idx].f!.ILC.daily_background_power
		+ integrator.f.f.vertices![idx].f!.ILC.kappa * integrator.f.f.vertices![idx].f!.ILC.mismatch_yesterday)
	end
	nothing
end
