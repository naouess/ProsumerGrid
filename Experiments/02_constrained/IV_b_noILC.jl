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
	u = [zeros(1000,1);1;zeros(1000,1)];
	fc = 1/6
	a = digitalfilter(Lowpass(fc), Butterworth(2))
	Q1 = filtfilt(a, u)
	Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);

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
	myPV = PV(ξ = t -> dd(t)[2], η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T_inv = 1/0.04), ILC = ILC_pars0, M_inv = 1/5)
	myLoad = Load(ξ = t -> dd(t)[1], η_load = t -> 1., LI = LI(kp = 110., ki = 0.004, T_inv = 1/0.045), ILC = ILC_pars1, M_inv = 1/4.8)
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
	ode_problem = ODEProblem(nd, i0, tspan)
end

@time sol = solve(ode_problem, Rodas4())#, maxiters=1e7)
# @show sol.destats

# Quick plots for validation
begin
	plot(sol, vars = syms_containing(nd, "ω"), legend = true)
	# savefig("$dir/plots/IV_real_noILC_omega.png")
end
plot(sol, vars = syms_containing(nd, "ϕ"), legend = true)
plot(sol, vars = syms_containing(nd, "integrated_LI"), legend = true)
plot(sol, vars = syms_containing(nd, "level"), legend = true)

## Extract hourly control powers
LI_exact, ILC_power = hourly_energy(sol, nd, num_days, N)
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
	    [level[:,1],LI_exact[:, 4]],
	    label = ["SoC" "LI"], #"ILC"],
		# legend=false,
		legendfontsize = 5,
		legend = :bottomright,
		lc=[:turquoise :orangered :black],
		# alpha=.5,
	    layout = (2, 1),
		xticks = (0:24:num_days*24*3600, string.(0:num_days)),
		ylims=[(0,1) (-0.4, 0.2) (-0.4, 0.2)],
		xaxis=("Days",font(8)),
		# yaxis=("Normed power",font(8))
		)
	# savefig("$dir/plots/IV_real_noILC_BatteryStuff.pdf")
end

begin
	psum = plot()
	plot!(0:3600:num_days*l_day, t -> dd(t)[1] - dd(t)[2], alpha=0.2, linewidth=3, linestyle=:dot) #.- dd(t)[2]
	plot!(1:3600:24*num_days*3600, (LI_exact[1:num_days*24,1] + LI_exact[1:num_days*24,2] + LI_exact[1:num_days*24,3] + LI_exact[1:num_days*24,4]),linewidth=3,
	linestyle=:dash,
	        xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
		     xtickfontsize=18,legend=false, legendfontsize=10,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm)
	# savefig("$dir/plots/IV_real_noILC_sum_exact.pdf")
end

integral = integrals(sol, nd, N, num_days, LI_exact, ILC_power)
sum_LI_n = integral[1]
sum_LI_p = integral[2]
sum_ILC_n = zeros(N)
sum_ILC_p = zeros(N)
sum_data_load = integral[5]
sum_data_infeed = integral[6]

LI_percentage = (sum(sum_LI_n) + sum(sum_LI_p)) ./ (sum(sum_LI_n) + sum(sum_LI_p) + sum(sum_ILC_n) + sum(sum_ILC_p))
LI_percentage_node = (sum_LI_n + sum_LI_p) ./ (sum_LI_n + sum_LI_p + sum_ILC_n + sum_ILC_p)

mismatch_data = (-(sum(sum_LI_n) + sum(sum_ILC_n) + sum_data_load) + (sum(sum_LI_p) + sum(sum_ILC_p) + sum_data_infeed)) / (sum_data_load +  sum_data_infeed)
mismatch_overall = (-(sum(sum_LI_n) + sum(sum_ILC_n) + sum_data_load) + (sum(sum_LI_p) + sum(sum_ILC_p) + sum_data_infeed)) / ((sum(sum_LI_n) + sum(sum_ILC_n) + sum_data_load) + (sum(sum_LI_p) + sum(sum_ILC_p) + sum_data_infeed))

indices_ω = idx_containing(nd, :ω)
KP = [nd.f.vertices![idx].f!.LI.kp for idx in 1:N]
maximum(maximum.(abs.(sol(t)[indices_ω[j]]) for j in 1:N for t in sol.t))
LI_exact_solt = zeros(num_days*24*3600, N)
data_solt = zeros(num_days*24*3600)
ILC_power_solt = zeros(num_days*24*3600)
for i in 1:num_days*24*3600
	for j in 1:N
		LI_exact_solt[i, j] = - KP[j] * sol(i)[indices_ω[j]] + sol(i)[indices_ω[j]+1]
	end
	data_solt[i] = dd(i)[1] - dd(i)[2]
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
	#fill=(0, 0, :lightblue),
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	linecolor=:black, linestyle=:dot,
	#linewidth=1
	)
	# savefig("$dir/plots/IV_real__noILC_mismatch_secondwise.pdf")
end


# Nodewise plots
begin
	node = 1
	p1 = plot()
	plot!(0:3600:num_days*l_day, t -> dd(t)[2], alpha=0.2, linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node],linewidth=3,
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10,  yaxis=("Normed power",font(10)),legend=false,  margin=5Plots.mm)
	ylims!(-1,1.)
	# savefig("$dir/plots/IV_real_noILC_node1.pdf")
end

begin
	node = 2
	p2 = plot()
	plot!(0:num_days*l_day, t -> dd(t)[1], alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node],linewidth=3,
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10,  yaxis=("Normed power",font(10)),legend=false,  margin=5Plots.mm)
    ylims!(-1,1)
	# savefig("$dir/plots/IV_real_noILC_node2.pdf")
end

begin
	node = 3
	p3 = plot()
	plot!(0:num_days*l_day, t -> 0, alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node],linewidth=3,
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10, yaxis=("Normed power",font(10)),legend=false, margin=5Plots.mm)
	ylims!(-1,1)
	# savefig("$dir/plots/IV_real_noILC_node3.pdf")
end

begin
	node = 4
	p4 = plot()
	plot!(0:num_days*l_day, t -> 0., alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, LI_exact[1:num_days*24, node], linewidth=3,
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=10,
	      xtickfontsize=10, legendfontsize=10, yaxis=("Normed power",font(10)),legend=false, margin=5Plots.mm)
	ylims!(-1,1)
	# savefig("$dir/plots/IV_real_noILC_node4.pdf")
end

begin
	l = @layout [a b; c d]
	plot_demand = plot(p1,p2,p3,p4, layout = l, xaxis=("Days", font(10)))
	# savefig("$dir/plots/IV_real_noILC_layout4Nodes.pdf")
end
