begin
	dir = @__DIR__
	include("$dir/src/ProsumerGrid.jl")
end

begin
	using ToeplitzMatrices
	using DSP
end

# General parameters
begin
	num_days = 10
	l_day = 24*3600
	l_hour = 3600
	N = 4
end

# Import demand functions
include("$dir/src/MinigridDemand.jl")


# Define ILC parameters for each node
begin
	vc = 1:N
	cover = Dict([v => [] for v in vc])
	Q = Matrix(1.0I, 24, 24)

	ILC_pars0 = ILC(kappa = 1/l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars1 = ILC(kappa = 1 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars2 = ILC(kappa = 1 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars3 = ILC(kappa = 1 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
end

# Define nodes, lines and power grid structure
begin
	dd = t -> demand_real(t)
	myPV = PV(ξ = t -> dd(t)[2], η_gen = t -> 1., LI = LI(kp = 100., ki = 0.05, T_inv = 1/0.04), ILC = ILC_pars0, M_inv = 1/5)
	myLoad = Load(ξ = t -> dd(t)[1], η_load = t -> 1., LI = LI( kp = 100., ki = 0.004, T_inv = 1/0.045), ILC = ILC_pars1, M_inv = 1/4.8)
	mySlack = Slack(η_gen = t -> 1., LI = LI(kp = 225., ki = 0.001, T_inv = 1/0.05), ILC = ILC_pars2, M_inv = 1/4.8)
	myBattery = Battery(η = η(gen=t->1., load=t->1.), LI = LI(kp = 200., ki = 0.001, T_inv = 1/0.043), ILC = ILC_pars3, M_inv = 1/4.8, C= 310 ./ 50)

	nodes = [constructor(myPV), constructor(myLoad), constructor(mySlack), constructor(myBattery)]

	myLine1 = PowerLine(from=1, to=2, K=6.)
	myLine2 = PowerLine(from=1, to=4, K=6.)
	myLine3 = PowerLine(from=3, to=2, K=6.)
	myLine4 = PowerLine(from=4, to=2, K=6.)
	myLine5 = PowerLine(from=3, to=4, K=6.)
	mylines = [myLine1, myLine2, myLine3, myLine4, myLine5]
	lines = [StaticEdge(f! = myLine, dim = 1) for myLine in mylines]

	g = SimpleGraph(N)
	add_edge!(g, 1, 2)
	add_edge!(g, 1, 4)
	add_edge!(g, 3, 2)
	add_edge!(g, 3, 4)
	add_edge!(g, 4, 2)

	nd = network_dynamics(nodes, lines, g)
end

begin
	tspan = (0., num_days*l_day)
	i0 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.7, 0.]

	# Define Periodic and Bounds watching Callbacks
	cb = CallbackSet(PeriodicCallback(HourlyUpdate, l_hour),
					 PeriodicCallback(DailyUpdate, l_day),
					 )

	# Define the ODE-Problem
	ode_problem = ODEProblem(nd, i0, tspan, callback = cb)
end

@time sol = solve(ode_problem, Rodas4())

# Quick plots for validation
begin
	plot(sol, vars = syms_containing(nd, "ω"), legend = true)
	# savefig("$dir/plots/VII_omega.png")
end
indices_ω = idx_containing(nd, :ω)
maximum(maximum.(abs.(sol(t)[indices_ω[j]]) for j in 1:N for t in sol.t))

plot(sol, vars = syms_containing(nd, "ϕ"), legend = true)
plot(sol, vars = syms_containing(nd, "integrated_LI"), legend = true)
plot(sol, vars = syms_containing(nd, "level"), legend = true)

## Extract hourly control powers
LI_exact, ILC_power = hourly_energy(sol, nd, num_days, N)
indices_ω = idx_containing(nd, :ω)

KP = [nd.f.vertices![idx].f!.LI.kp for idx in 1:N]
for i in 1:24*num_days+1
	for j in 1:2
		LI_exact[i, j] = - KP[j] * sol((i-1)*3600)[indices_ω[j]] # + sol((i-1)*3600)[indices_ω[j]+1]
	end
end
indices = idx_containing(nd, :integrated_LI)
hourly = zeros(24 * num_days + 1, N)
for j = 1:N
	for i = 1:24*num_days+1
		hourly[i, j] = sol((i-1)*3600)[indices[j]]
	end
end

# Plot battery stuff
begin
	indices1 = idx_containing(nd, :level)
	vars = syms_containing(nd, :level)

	level = zeros(24 * num_days + 1, length(vars))
	for j = 1:length(vars)
		for i = 1:24*num_days+1
			level[i, j] = sol((i-1)*3600)[indices1[j]]
		end
	end
end

begin
	plot(
	    [level[:,1],LI_exact[:, 4], vcat(ILC_power[1:10,:, 4]'...)],
	    label = ["SoC" "LI" "ILC"],
		legend=false,
		legendfontsize = 5,
		#legend = :bottomright,
		lc=[:turquoise :orangered :black],
		# alpha=.5,
	    layout = (3, 1),
		xticks = (0:24:num_days*24*3600, string.(0:num_days)),
		ylims=[(0,1) (-0.4, 0.2) (-0.4, 0.2)],
		xaxis=("Days",font(8)),
		title = ["SoC" "LI" "ILC"], titleloc = :left, titlefont = font(8),
		# margin=5Plots.mm
		# yaxis=("Normed power",font(8))
		)
	# savefig("$dir/plots/VII_BatteryStuff.pdf")
end
batt = [battery(i) for i in 5*3600*24:3600:24*3600*9] ./ 310
plot(level[5*24:9*24-1,1],
	 xticks = (5:24:num_days*24*3600, string.(5:num_days)),
	 ylims=(0.2, 1),
	 xaxis=("Days",font(8)),
	 legend=false,
	 label = "Simulated",
	 yaxis=("State of charge in %",font(8)),
	 margin=5Plots.mm,
	 linewidth = 1.5,
)
plot!(batt, linewidth=1.5, linestyle=:dash, label="Real")
savefig("$dir/plots/VII_battery")

plot(0:num_days*l_day, t -> dd(t)[1] - dd(t)[2])

begin
	psum = plot()
	ILC_power_hourly_mean_sum = vcat(ILC_power[:,:,1]'...) .+ vcat(ILC_power[:,:,2]'...) .+ vcat(ILC_power[:,:,3]'...) .+ vcat(ILC_power[:,:,4]'...)
	plot!(0:3600:num_days*l_day, t -> dd(t)[1] - dd(t)[2], alpha=0.7, linewidth=3, linestyle=:dot) #.- dd(t)[2]
	plot!(1:3600:24*num_days*3600, (LI_exact[1:num_days*24,1] + LI_exact[1:num_days*24,2] + LI_exact[1:num_days*24,3] + LI_exact[1:num_days*24,4]),linewidth=3, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	               xtickfontsize=18,legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)),lc =:black, margin=5Plots.mm)
	savefig("$dir/plots/VII_sum_exact.pdf")
end

begin
	psum = plot()
	ILC_power_hourly_mean_sum = vcat(ILC_power[:,:,1]'...) .+ vcat(ILC_power[:,:,2]'...) .+ vcat(ILC_power[:,:,3]'...) .+ vcat(ILC_power[:,:,4]'...)
	plot!(0:3600:num_days*l_day, t -> dd(t)[1] - dd(t)[2], alpha=0.2, linewidth=3, linestyle=:dot) #.- dd(t)[2]
	plot!(1:3600:24*num_days*3600, (hourly[1:num_days*24,1] + hourly[1:num_days*24,2] + hourly[1:num_days*24,3] + hourly[1:num_days*24,4]) ./ 3600,linewidth=3)#, linestyle=:auto)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	               xtickfontsize=18,legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(14)),  yaxis=("Normed power",font(14)),lc =:black, margin=5Plots.mm)
	# savefig("$dir/plots/VI_allILC_sum_hourly.pdf")
end

begin
	control =  ILC_power_hourly_mean_sum[1:num_days*24] + (LI_exact[1:num_days*24,1] + LI_exact[1:num_days*24,2] + LI_exact[1:num_days*24,3] + LI_exact[1:num_days*24,4]) #./3600
	data = [- dd(i)[1] + dd(i)[2] for i in 0:3600:24*num_days*3600-1] #
	mismatch = data + control
	maximum(abs.(mismatch))
	# sum(abs.(mismatch)) / sum(abs.(data))
end


begin
	# plot exact power balance
	p_balance = plot()
	plot!(1:24*num_days*3600, t -> dd(t)[1] - dd(t)[2] , alpha=.9, fill=(0, 0.7, :lightblue), linecolor=:lightblue, linestyle=:dot,
	xticks = (0:24*3600:num_days*24*3600, string.(0:num_days)), label="Microgrid data",
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(14)),  yaxis=("Normed power",font(14)), margin=5Plots.mm)
	plot!(1:3600:24*num_days*3600, -mismatch, linewidth = 0, legend = false, linecolor=:black, alpha=1, fill=(0,.5,:black), label="Power balance difference")
	# savefig("$dir/plots/VI_allILC_mismatch_hourly.pdf")
end

integral = integrals(sol, nd, N, num_days, LI_exact, ILC_power)
sum_LI_n = integral[1]
sum_LI_p = integral[2]
sum_ILC_n = integral[3]
sum_ILC_p = integral[4]
sum_data_load = integral[5]
sum_data_infeed = integral[6]

LI_percentage = (sum(sum_LI_n) + sum(sum_LI_p)) ./ (sum(sum_LI_n) + sum(sum_LI_p) + sum(sum_ILC_n) + sum(sum_ILC_p))
LI_percentage_node = (sum_LI_n + sum_LI_p) ./ (sum_LI_n + sum_LI_p + sum_ILC_n + sum_ILC_p)

mismatch_data = (-(sum(sum_LI_n) + sum(sum_ILC_n) + sum_data_load) + (sum(sum_LI_p) + sum(sum_ILC_p) + sum_data_infeed)) / (sum_data_load +  sum_data_infeed)
mismatch_overall = (-(sum(sum_LI_n) + sum(sum_ILC_n) + sum_data_load) + (sum(sum_LI_p) + sum(sum_ILC_p) + sum_data_infeed)) / ((sum(sum_LI_n) + sum(sum_ILC_n) + sum_data_load) + (sum(sum_LI_p) + sum(sum_ILC_p) + sum_data_infeed))

indices_ω = idx_containing(nd, :ω)
KP = [nd.f.vertices![idx].f!.LI.kp for idx in 1:N]

LI_exact_solt = zeros(num_days*24*3600, N)
data_solt = zeros(num_days*24*3600)
ILC_power_solt = zeros(num_days*24*3600)
for i in 1:num_days*24*3600
	for j in 1:N
		LI_exact_solt[i, j] = - KP[j] * sol(i)[indices_ω[j]] + sol(i)[indices_ω[j]+1]
	end
	data_solt[i] = dd(i)[1] - dd(i)[2]
	hour = Int(div(i, 3600) + 1)
	ILC_power_solt[i] = ILC_power_hourly_mean_sum[hour]
end

control_solt = ILC_power_solt + LI_exact_solt[:, 1] + LI_exact_solt[:, 2] + LI_exact_solt[:, 3] + LI_exact_solt[:, 4]
mismatch_solt = - data_solt + control_solt

begin
	plot(data_solt, labels="Data", linecolor=:lightblue,
	fill=(0, 0, :lightblue),
	# alpha=1,
	#linestyle=:dot,
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm,
	# ribbon=mismatch_solt
	)
	plot!(control_solt,
	#mismatch_solt,
	#fill=(0, 0, :lightblue),
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	linecolor=:black, linestyle=:dot,
	#linewidth=1
	)
	# savefig("$dir/plots/VII_mismatch_secondwise.png")
end


# Nodewise plots
using Plots
begin
	node = 1
	p1 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:3600:num_days*l_day, t -> dd(t)[2], alpha=0.7, linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node],linewidth=3,
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10, yaxis=("Normed power",font(10)),legend=false, margin=5Plots.mm)
	ylims!(-1,1.)
	savefig("$dir/plots/VII_node1.png")
end

begin
	node = 2
	p2 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> dd(t)[1], alpha=0.7,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node],linewidth=3,
	# plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24],
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10, #linewidth=3, lc =:black,
	      xtickfontsize=10, legendfontsize=10, yaxis=("Normed power",font(10)),legend=false, margin=5Plots.mm)
    ylims!(-1,1)
	savefig("$dir/plots/VII_node2.png")
end

begin
	node = 3
	p3 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> 0, alpha=0.7,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node],linewidth=3)
	plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24], xticks = (1:3600*24:num_days*24*3600+1, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10, linewidth=3, yaxis=("Normed power",font(10)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-1,1)
	savefig("$dir/plots/VII_node3.png")
end

begin
	node = 4
	p4 = plot()
	ILC_power_hourly_mean_node = vcat(ILC_power[:,:,node]'...)
	plot!(0:num_days*l_day, t -> 0., alpha=0.7,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node], linewidth=3)
	plot!(1:3600:num_days*24*3600, ILC_power_hourly_mean_node[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10, linewidth=3, yaxis=("Normed power",font(10)),legend=false, lc =:black, margin=5Plots.mm)
	ylims!(-1,1)
	savefig("$dir/plots/VII_node4.png")
end

node = 3
ILC_power_node = zeros(length(sol.t))
LI_power_node = zeros(length(sol.t))
control_node = zeros(length(sol.t))
control_integral_node = 0

indices_ω = idx_containing(nd, :ω)
indices_χ = idx_containing(nd, :χ)
KP = [nd.f.vertices![idx].f!.LI.kp for idx in 1:N]

for (i,t) in enumerate(sol.t)
	if i == 1
		delta_t = (sol.t[i+1] - sol.t[i])
	elseif i == length(sol.t)
		delta_t = (sol.t[i] - sol.t[i-1])
	else
		delta_t = (sol.t[i+1] - sol.t[i-1])/2
	end

	ILC_node = vcat(ILC_power[:,:,node]'...)
	hour = Int(div(t, 3600) + 1)

	ILC_power_node[i] = ILC_node[hour]
	LI_power_node[i] = - KP[node] * sol(t)[indices_ω[node]] + sol(t)[indices_χ[node]]

	control_node[i] = LI_power_node[i] + ILC_power_node[i]

	global control_integral_node += abs(control_node[i]) * delta_t
	# global control_solt_integral += abs(control_solt[i])
end

plot(control_node, 	xticks=(0:3600:l_day*num_days, string.(0:num_days)))
control_integral_node / sum(abs.(control_solt))
(sum_LI_p[node] - sum_LI_n[node])/ sum(abs.(control_solt))


observables = integrals_from_to(sol, nd, N, num_days, LI_exact, ILC_power)
C_plus_diesel = (observables[2][3] + observables[4][3]) / (observables[5] + observables[6])
C_minus_diesel = (observables[1][3] + observables[3][3]) / (observables[5] + observables[6])
C_diesel_real = (sum([diesel(i) ./ 50. for i in 5*24:24*num_days*3600])) / (observables[5] + observables[6])
