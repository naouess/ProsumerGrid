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

function watch_storage_1(network)

	indices = idx_containing(network, :level)
	# storage_nodes = syms_containing(network, :level)
	# indices_control = idx_containing(network, :integrated_LI)

    function condition(out, u, t, integrator)
		out[2] = u[indices[1]] - 1.
		out[1] = u[indices[1]]
    end

    function affect!(integrator, idx)
		@info "bounds touched at t=$(integrator.t)."
		if idx == 1
			integrator.u[indices[1]] = 0.
		end
		if idx == 2
			integrator.u[indices[1]] = 1.
		end
		# integrator.f.f.vertices![4].f!.ILC.current_background_power = - 0.98 * integrator.f.f.vertices![4].f!.ILC.current_background_power
		# integrator.u[indices[1] - 1] = - 0.98 * integrator.u[indices[1] - 1]
    end
    return VectorContinuousCallback(condition, affect!, 2)
end

#
# @doc """
#     watch_storage_discharge(network)
# Watching the discharge lower bound. Currently set to 0. Could be adjusted to a more realistic value depending on storage type.
# """
# function watch_storage_discharge(network)
#
# 	indices = idx_containing(network, :level)
# 	storage_nodes = syms_containing(network, :level)
# 	# indices_control = idx_containing(network, :integrated_LI)
#
#     function condition(u, t, integrator)
# 		# true
# 		return u[indices[1]] <= 0.
#     end
#
#     function affect!(integrator)
# 		# if integrator.u[indices[1]] <= 0.
# 		# @info "Discharging bounds touched at t=$(integrator.t)."
# 		integrator.f.f.vertices![4].f!.ILC.current_background_power = 0.
# 		# integrator.u[indices[1]] = 0.
# 		# println(integrator(integrator.t, Val{1}).x)
# 		# integrator(integrator.t, Val{1}).x[1][1] = 0. #indices[1]-1
# 		# integrator.f.f.vertices![4].f!.ILC.current_background_power = 0.
# 		# end
#     end
#     return DiscreteCallback(condition, affect!, save_positions = (false,false))
# end
#
# @doc """
#     watch_storage_charge(network)
# Watching the discharge lower bound. Currently set to 0. Could be adjusted to a more realistic value depending on storage type.
# """
# function watch_storage_charge(network)
# 	indices = idx_containing(network, :level)
# 	storage_nodes = syms_containing(network, :level)
#
#     function condition(u, t, integrator)
# 		# true
# 		return u[indices[1]] >= 1.
#     end
#
#     function affect!(integrator)
# 		# if integrator.u[indices[1]] >= 1.
# 		# @info "Charging bounds touched at t=$(integrator.t)."
# 		integrator.f.f.vertices![4].f!.ILC.current_background_power = 0.
# 		# integrator.u[indices[1]] = 1.
# 		# integrator(integrator.t, Val{1}) = 0.
# 		# integrator.f.f.vertices![4].f!.ILC.current_background_power = 0.
# 		# end
#     end
#     return DiscreteCallback(condition, affect!, save_positions = (false,false))
# end

## SavingCallback
# function plot_flow(saved_edge_data, network)
#     t = saved_edge_data.t
#     vals = saved_edge_data.saveval
#
#     data = zeros(length(t), length(vals[1]))
#     for i in 1:length(vals)
#         data[i, :] = abs.(vals[i])
#     end
#     plot(t, data, size=(800,600))
# end
#
# # define callback to save edge values (for plotting)
# function save_edges(u, t, integrator)
#     graph_data = integrator.f.f(u, integrator.p, integrator.t, GetGD)
#     return copy(graph_data.v_array)
# end
