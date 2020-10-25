begin
	dir = @__DIR__
	include("$dir/src/ProsumerGrid_withoutConstraints.jl")
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
	real_topology = true
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

	ILC_pars0 = ILC(kappa = 0.75 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars1 = ILC(kappa = 0.75 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars2 = ILC(kappa = 0.75 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
	ILC_pars3 = ILC(kappa = 0.75 /l_hour, mismatch_yesterday = zeros(24), daily_background_power = zeros(24),
	current_background_power = 0., ilc_nodes = vc, ilc_cover = cover, Q = Q)
end

# Define nodes, lines and power grid structure
begin
	dd = t -> demand_real(t)
	myPV = PV(ξ = t -> dd(t)[2], η_gen = t -> 1., LI = LI(kp = 400., ki = 0.05, T_inv = 1/0.04), ILC = ILC_pars0, M_inv = 1/5)
	myLoad = Load(ξ = t -> dd(t)[1], η_load = t -> 1., LI = LI( kp = 110., ki = 0.004, T_inv = 1/0.045), ILC = ILC_pars1, M_inv = 1/4.8)
	mySlack = Slack(η_gen = t -> 1., LI = LI(kp = 100., ki = 0.05, T_inv = 1/0.047), ILC = ILC_pars2, M_inv = 1/4.1)
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
	i0 = zeros(21)

	# Define Periodic and Bounds watching Callbacks
	cb = CallbackSet(PeriodicCallback(HourlyUpdate, l_hour),
					 PeriodicCallback(DailyUpdate, l_day))

	# Define the ODE-Problem
	ode_problem = ODEProblem(nd, i0, tspan)#, callback = cb)
end

@time sol = solve(ode_problem, Rodas4())

begin
	plot(sol, vars = syms_containing(nd, "ω"), legend = true)
	# savefig("$dir/plots/III_real_omega.png")
end

indices_ω = idx_containing(nd, :ω)
maximum(maximum.(abs.(sol(t)[indices_ω[j]]) for j in 1:N for t in sol.t))

# Extract hourly control powers
LI_exact, ILC_power = hourly_energy(sol, nd, num_days, N)
indices = idx_containing(nd, :integrated_LI)
ILC_power_hourly_mean_sum = zeros(24 * num_days + 1)

begin
	psum = plot()
	plot!(0:3600:num_days*l_day, t -> dd(t)[1] - dd(t)[2], alpha=0.2, linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, (LI_exact[1:num_days*24,1] + LI_exact[1:num_days*24,2] + LI_exact[1:num_days*24,3] + LI_exact[1:num_days*24,4]),
	linewidth=3, linestyle=:dash,
	xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), ytickfontsize=14,
	xtickfontsize=18,legend=false, legendfontsize=10,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)),margin=5Plots.mm)
	# savefig("$dir/plots/III_real_sum_exact.pdf")
end

begin
	control =  ILC_power_hourly_mean_sum[1:num_days*24] + (LI_exact[1:num_days*24,1] + LI_exact[1:num_days*24,2] + LI_exact[1:num_days*24,3] + LI_exact[1:num_days*24,4]) #./3600
	data = [- dd(i)[1] + dd(i)[2] for i in 0:3600:24*num_days*3600-1] #
	mismatch = data + control
	maximum(abs.(mismatch))
end

begin
	# plot exact power balance
	p_balance = plot()
	plot!(1:24*num_days*3600, t -> dd(t)[1] - dd(t)[2] , alpha=.9, fill=(0, 0.7, :lightblue), linecolor=:lightblue, linestyle=:dot,
	xticks = (0:24*3600:num_days*24*3600, string.(0:num_days)), label="Microgrid data",
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(14)),  yaxis=("Normed power",font(14)), margin=5Plots.mm)
	plot!(1:3600:24*num_days*3600, -mismatch, linewidth = 0, legend = false, linecolor=:black, alpha=1, fill=(0,.5,:black), label="Power balance difference")
	# savefig("$dir/plots/III_real_mismatch_hourly.pdf")
end
ILC_power = zeros(num_days+2, 24, N)
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
	# mismatch_solt,
	#fill=(0, 0, :lightblue),
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	linecolor=:black,
	linestyle=:dot,
	linewidth=1
	)
	# savefig("$dir/plots/III_real_mismatch_secondwise.pdf")
end
