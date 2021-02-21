## Net flow, from node j to all neighbouring nodes k
@inline function total_flow(e_s, e_d)
	net_flow = 0.
    @inbounds for e in e_s
        net_flow -= e[1]
    end
    @inbounds for e in e_d
        net_flow += e[1]
    end
	net_flow
end

## TODO:
#  Define useful extracting helpers
@inline function hourly_energy(sol, nd, num_days, N)
	indices = idx_containing(nd, :integrated_LI)
	hourly_energy = zeros(24 * num_days + 1, N)
	for j = 1:N
		for i = 1:24*num_days+1
			hourly_energy[i, j] = sol((i-1)*3600)[indices[j]]
		end
	end
	ILC_power = zeros(num_days+2, 24, N)
	# Values for ILC_power for day 1:
	kappa = nd.f.vertices![1].f!.ILC.kappa
	Q_filter = nd.f.vertices![1].f!.ILC.Q

	for j = 1:N
		ILC_power[2,:,j] = Q_filter * (zeros(24,1) +  kappa * hourly_energy[1:24, j])
	end
	# Values for ILC_power from day 2 until end:
	for i= 2:num_days
		for j = 1:N
			ILC_power[i+1,:,j] .= Q_filter * (ILC_power[i,:,j] .+  kappa*hourly_energy[(i-1)*24+1:i*24,j])
		end
	end

	indices_ω = idx_containing(nd, :ω)
	KP = [nd.f.vertices![idx].f!.LI.kp for idx in 1:N]
	LI_exact = zeros(24 * num_days + 1, N)
	for i in 1:24*num_days+1
		for j in 1:N
			LI_exact[i, j] = - KP[j] * sol((i-1)*3600)[indices_ω[j]] + sol((i-1)*3600)[indices_ω[j]+1]
		end
	end
	return LI_exact, ILC_power
end

## Observables
@inline function integrals(sol, nd, N, num_days, LI_exact, ILC_power)
	sum_LI_n = zeros(N)
	sum_LI_p = zeros(N)
	sum_ILC_n = zeros(N)
	sum_ILC_p = zeros(N)

	# sum_mismatch = 0.
	sum_data_load = 0.
	sum_data_infeed = 0.
	ILC_power_hourly_mean_sum = zeros((num_days+2)*24)

	indices_ω = idx_containing(nd, :ω)
	indices_χ = idx_containing(nd, :χ)
	KP = [nd.f.vertices![idx].f!.LI.kp for idx in 1:N]

	for (i,t) in enumerate(sol.t) #sol.t #1:24*num_days+1
		if i == 1
			delta_t = (sol.t[i+1] - sol.t[i])
		elseif i == length(sol.t)
			delta_t = (sol.t[i] - sol.t[i-1])
		else
			delta_t = (sol.t[i+1] - sol.t[i-1])/2
		end
		for j in 1:N
			if sign(- KP[j] * sol(t)[indices_ω[j]] + sol(t)[indices_χ[j]]) == -1
				sum_LI_n[j] += abs(- KP[j] * sol(t)[indices_ω[j]] + sol(t)[indices_χ[j]]) * delta_t
			elseif sign(- KP[j] * sol(t)[indices_ω[j]] + sol(t)[indices_χ[j]]) == 1
				sum_LI_p[j] += abs(- KP[j] * sol(t)[indices_ω[j]] + sol(t)[indices_χ[j]]) * delta_t
			end
		end
		sum_data_load = sum_data_load + abs(dd(t)[1]) * delta_t
		sum_data_infeed = sum_data_infeed  + abs(dd(t)[2]) * delta_t
		# LI_demand = sum_LI[j] - sum_data TODO hier definieren
	end

	for j in 1:N
 		ILC = vcat(ILC_power[:,:,j]'...)
		ILC_power_hourly_mean_sum = ILC_power_hourly_mean_sum .+ ILC
		for i in 1:24*num_days
			if sign(ILC[i]) == -1
				sum_ILC_n[j] += abs(ILC[i]) * 3600
			elseif sign(ILC[i]) == 1
				sum_ILC_p[j] += abs(ILC[i]) * 3600
			end
		end
	end
	return sum_LI_n, sum_LI_p, sum_ILC_n, sum_ILC_p, sum_data_load, sum_data_infeed #, sum_1
end

## TODO:
#  Define useful plotting helpers

# using LaTeXStrings
# using Plots, GraphPlot
# using LinearAlgebra
# using Statistics

# @inline function plot_graph()
# 	gplot()
# end
#
# @inline function plot_frequency(node, time)
# 	plot(sol, vars = syms_containing(nd, "ω"), legend = true)
# end
#
# @inline function plot_energy()
# 	plot(sol, vars = syms_containing(nd, "integrated_LI"), legend = true)
# end
#
# @inline function all_plots()
#
# 	# TODO here add option (savefig or not)
#
# 	# NODE WISE second-wisenode = 1
# 	begin
# 		node = 1
# 		p1 = plot()
# 		ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
# 		plot!(0:num_days*l_day, t -> dd(t)[node], alpha=0.2, label = latexstring("P^d_$node"),linewidth=3, linestyle=:dot)
# 		plot!(1:3600:24*num_days*3600,hourly_energy[1:num_days*24, node]./3600, label=latexstring("y_$node^{c,h}"),linewidth=3)
# 		plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_node[1:num_days*24], label=latexstring("\$u_$node^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
# 		      xtickfontsize=14, legendfontsize=10, linewidth=3, yaxis=("normed power",font(14)),legend=false, lc =:black, margin=5Plots.mm)
# 		ylims!(-0.7,1.5)
# 		title!(L"j = 1")
# 		savefig("$dir/plots/demand_seconds_node_$(node)_hetero.png")
# 	end
#
# 	psum = plot()
# 	ILC_power_hourly_mean_sum = vcat(ILC_power[:,:,1]'...) .+ vcat(ILC_power[:,:,2]'...) .+ vcat(ILC_power[:,:,3]'...) .+ vcat(ILC_power[:,:,4]'...)
# 	plot!(0:num_days*l_day, t -> (dd(t)[1] .+ dd(t)[2] .+ dd(t)[3] .+ dd(t)[4]), alpha=0.2, label = latexstring("\$P^d_j\$"),linewidth=3, linestyle=:dot)
# 	plot!(1:3600:24*num_days*3600,(hourly_energy[1:num_days*24,1] + hourly_energy[1:num_days*24,2] + hourly_energy[1:num_days*24,3] + hourly_energy[1:num_days*24,4])./3600, label=latexstring("y_j^{c,h}"),linewidth=3, linestyle=:dash)
# 	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], label=latexstring("\$u_j^{ILC}\$"), xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
# 	               xtickfontsize=18,legend=false, legendfontsize=10, linewidth=3,xaxis=("days [c]",font(14)),  yaxis=("normed power",font(14)),lc =:black, margin=5Plots.mm)
# 	# ylims!(-0.7,1.5)
# 	# title!("Initial convergence")
# 	savefig(psum,"$dir/plots/demand_seconds_sum_hetero.png")
# end
#
# plot(sol, vars = syms_containing(nd, "ϕ"), legend = true)
# plot(sol, vars = syms_containing(nd, "integrated_LI"), legend = true)
# plot(sol, vars = syms_containing(nd, "level"), legend = true)
