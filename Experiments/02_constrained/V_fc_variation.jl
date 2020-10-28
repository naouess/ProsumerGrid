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
	f_c_lst = 1 ./ 1 #(1:1:24)
	f_c = 1/2
	num_monte = length(f_c_lst)

	a = digitalfilter(Lowpass(f_c), Butterworth(2))
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
	add_edge!(g, 2, 4)

	nd = network_dynamics(nodes, lines, g)
end

begin
	tspan = (0., num_days*l_day)
	i0 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.7, 0.]

	# Define Periodic and Bounds watching Callbacks
	cb = CallbackSet(PeriodicCallback(HourlyUpdate, l_hour),
					 PeriodicCallback(DailyUpdate, l_day))

	# Define the ODE-Problem
	ode_problem = ODEProblem(nd, i0, tspan, callback = cb)
end

batch_size = num_monte
get_run(i, batch_size) = mod(i, batch_size)==0 ? batch_size : mod(i, batch_size)

function prob_func_fc(prob, i, repeat, batch_size, f_c_lst, num_days, N, cb)
	println("sim ", i)
	run = get_run(i, batch_size)

	if f_c_lst[run] == 1.
		Q = Matrix(1.0I, 24, 24)
	else
		u = [zeros(1000,1);1;zeros(1000,1)];
		a = digitalfilter(Lowpass(f_c_lst[run]), Butterworth(2))
		Q1 = filtfilt(a, u)
		Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);
	end
	for idx in 1:N
		prob.f.f.vertices![idx].f!.ILC.daily_background_power = zeros(24)
		prob.f.f.vertices![idx].f!.ILC.current_background_power = 0
		prob.f.f.vertices![idx].f!.ILC.mismatch_yesterday = zeros(24)
		prob.f.f.vertices![idx].f!.ILC.Q = Q
	end

	ODEProblem(prob.f, prob.u0, prob.tspan, callback=cb)
	prob
end

# Observer function: making sense of solution: defining observables and key indicators
function observer_fc(sol, i, nd, num_days, N)
	LI_exact, ILC_power = hourly_energy(sol, sol.prob.f, num_days, N)
	ILC_power_hourly_mean_sum = vcat(ILC_power[:,:,1]'...) .+ vcat(ILC_power[:,:,2]'...) .+ vcat(ILC_power[:,:,3]'...) .+ vcat(ILC_power[:,:,4]'...)

	indices_ω = idx_containing(nd, :ω)
	omega_max = maximum(maximum.(abs.(sol(t)[indices_ω[j]]) for j in 1:N for t in sol.t))
	all_integrals = integrals(sol, sol.prob.f, N, num_days, LI_exact, ILC_power)

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

	((LI_exact, ILC_power, all_integrals, sol.prob.f.f.vertices![1].f!.ILC.Q, sol.t, omega_max, LI_exact_solt, ILC_power_solt, data_solt), false)
end

monte_prob = EnsembleProblem(
	ode_problem,
	output_func = (sol, i) -> observer_fc(sol, i, nd, num_days, N),
	prob_func = (prob,i,repeat) -> prob_func_fc(prob, i, repeat, batch_size, f_c_lst, num_days, N, cb),
	u_init = [])

@time res = solve(monte_prob,
			      Rodas4P(),
				  trajectories=num_monte,
				  )


Q_filter = [p[4] for p in res.u]

# to validate that Q filter was correctly set
for (i, f_c) in enumerate(f_c_lst)
    if f_c == 1
        Q_valid = Matrix(1.0I, 24, 24)
    else
        a = digitalfilter(Lowpass(f_c), Butterworth(2))
        Q1 = filtfilt(a, u)
        Q = Toeplitz(Q1[1001:1001+24-1], Q1[1001:1001+24-1]);
        Q_valid = Q
    end
    println(i)
    @assert Q_valid == Q_filter[i]
end

LI_integral_n = [p[3][1] for p in res.u]
LI_integral_p = [p[3][2] for p in res.u]

ILC_integral_n = [p[3][3] for p in res.u]
ILC_integral_p = [p[3][4] for p in res.u]

data_integral_load = [p[3][5] for p in res.u]
data_integral_infeed = [p[3][6] for p in res.u]

t_res = [p[5] for p in res.u]

LI_exact = [p[1] for p in res.u]
ILC_power = [p[2] for p in res.u]
LI_exact_solt = [p[7] for p in res.u]
ILC_power_solt = [p[8] for p in res.u]
data_solt = [p[9] for p in res.u]

LI_sum_n = zeros(num_monte)
ILC_sum_n = zeros(num_monte)
LI_sum_p = zeros(num_monte)
ILC_sum_p = zeros(num_monte)
percentage_LI_sum = zeros(num_monte)
omega_max = [p[6] for p in res.u]

for k in 1:num_monte
	for j in 1:N
		LI_sum_n[k] += LI_integral_n[k][j]
		LI_sum_p[k] += LI_integral_p[k][j]
		ILC_sum_n[k] += ILC_integral_n[k][j]
		ILC_sum_p[k] += ILC_integral_p[k][j]
	end
	percentage_LI_sum[k] = (LI_sum_n[k] + LI_sum_p[k])/ (LI_sum_n[k] + LI_sum_p[k] + ILC_sum_p[k] + ILC_sum_n[k])
end
percentage_LI_sum
mismatch_overall = ((LI_sum_p + ILC_sum_p + data_integral_infeed) - (LI_sum_n + ILC_sum_n + data_integral_load)) ./ ((LI_sum_p + ILC_sum_p + data_integral_infeed) + (LI_sum_n + ILC_sum_n + data_integral_load))
mismatch_data = ((LI_sum_p + ILC_sum_p + data_integral_infeed) - (LI_sum_n + ILC_sum_n + data_integral_load)) ./ (data_integral_infeed + data_integral_load)

using Plots
using LaTeXStrings
markershapes= [:circle];
markercolors= [:orange];
plot(mismatch_overall*100,
	 xticks = (1:1:24),
	 xaxis = ("Simulation", font(10)),
	 yaxis = ("Power balance difference in %", font(10)),
	 # ylabel = L"\Delta_{%, overall}",
	 legend = false,
	 margin=5Plots.mm,
	 linestyle = :solid,
	 ylims = (-2, -0.5),
	 yticks = (-2:0.5:-0.5),
	 linewidth = 1,
	 color = markercolors,
	 shape = markershapes,
	 )
# savefig("$dir/plots/V_Mismatch.png")

plot(percentage_LI_sum*100,
	 xticks = (1:1:24),
	 xaxis = ("Simulation", font(10)),
	 yaxis = ("Share of LI power in %", font(10)),
	 # ylabel = L"\Delta_{%, overall}",
	 legend = false,
	 margin=5Plots.mm,
	 linestyle = :solid,
	 ylims = (40, 92),
	 # yticks = (-2:0.5:-0.5),
	 linewidth = 1,
	 color = markercolors,
	 shape = markershapes,
	 )
# savefig("$dir/plots/V_LI.png")

s = 2
begin
	psum_2 = plot()
	ILC_power_hourly_mean_sum =  vcat(ILC_power[s][:,:,1]'...) .+ vcat(ILC_power[s][:,:,2]'...) .+ vcat(ILC_power[s][:,:,3]'...) .+ vcat(ILC_power[s][:,:,4]'...)
	plot!(0:num_days*l_day, t -> (dd(t)[1] .- dd(t)[2]), alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, (LI_exact[s][1:num_days*24,1] + LI_exact[s][1:num_days*24,2] + LI_exact[s][1:num_days*24,3] + LI_exact[s][1:num_days*24,4]),linewidth=3, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), #ytickfontsize=14, xtickfontsize=18,
	               legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), lc =:black, margin=5Plots.mm)
	# savefig("$dir/plots/V_Qfilter2_sum_exact.pdf")
end

control_solt = ILC_power_solt[s] + LI_exact_solt[s][:, 1] + LI_exact_solt[s][:, 2] + LI_exact_solt[s][:, 3] + LI_exact_solt[s][:, 4]
mismatch_solt = - data_solt[s] + control_solt

begin
	plot(data_solt[s], labels="Data", linecolor=:lightblue,
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
	# savefig("$dir/plots/V_Qfilter2_ILC_mismatch_secondwise.pdf")
end

s = 24
begin
	psum_24 = plot()
	ILC_power_hourly_mean_sum =  vcat(ILC_power[s][:,:,1]'...) .+ vcat(ILC_power[s][:,:,2]'...) .+ vcat(ILC_power[s][:,:,3]'...) .+ vcat(ILC_power[s][:,:,4]'...)
	plot!(0:num_days*l_day, t -> (dd(t)[1] .- dd(t)[2]), alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, (LI_exact[s][1:num_days*24,1] + LI_exact[s][1:num_days*24,2] + LI_exact[s][1:num_days*24,3] + LI_exact[s][1:num_days*24,4]),linewidth=3, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), #ytickfontsize=14, xtickfontsize=18,
	               legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), lc =:black, margin=5Plots.mm)
	# savefig("$dir/plots/V_Qfilter24_sum_exact.pdf")
end

control_solt = ILC_power_solt[s] + LI_exact_solt[s][:, 1] + LI_exact_solt[s][:, 2] + LI_exact_solt[s][:, 3] + LI_exact_solt[s][:, 4]
mismatch_solt = - data_solt[s] + control_solt

begin
	plot(data_solt[s], labels="Data", linecolor=:lightblue,
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
	# savefig("$dir/plots/V_Qfilter24_mismatch_secondwise.pdf")
end

s = 6
begin
	psum_6 = plot()
	ILC_power_hourly_mean_sum =  vcat(ILC_power[s][:,:,1]'...) .+ vcat(ILC_power[s][:,:,2]'...) .+ vcat(ILC_power[s][:,:,3]'...) .+ vcat(ILC_power[s][:,:,4]'...)
	plot!(0:num_days*l_day, t -> (dd(t)[1] .- dd(t)[2]), alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, (LI_exact[s][1:num_days*24,1] + LI_exact[s][1:num_days*24,2] + LI_exact[s][1:num_days*24,3] + LI_exact[s][1:num_days*24,4]),linewidth=3)#, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), #ytickfontsize=14, xtickfontsize=18,
	               legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), lc =:black, margin=5Plots.mm)
	# savefig("$dir/plots/V_Qfilter6_sum_exact.pdf")
end

control_solt = ILC_power_solt[s] + LI_exact_solt[s][:, 1] + LI_exact_solt[s][:, 2] + LI_exact_solt[s][:, 3] + LI_exact_solt[s][:, 4]
mismatch_solt = - data_solt[s] + control_solt

begin
	plot(data_solt[s], labels="Data", linecolor=:lightblue,
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
	# savefig("$dir/plots/V_Qfilter6_mismatch_secondwise.pdf")
end

s = 12
begin
	psum_12 = plot()
	ILC_power_hourly_mean_sum =  vcat(ILC_power[s][:,:,1]'...) .+ vcat(ILC_power[s][:,:,2]'...) .+ vcat(ILC_power[s][:,:,3]'...) .+ vcat(ILC_power[s][:,:,4]'...)
	plot!(0:num_days*l_day, t -> (dd(t)[1] .- dd(t)[2]), alpha=0.2,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, (LI_exact[s][1:num_days*24,1] + LI_exact[s][1:num_days*24,2] + LI_exact[s][1:num_days*24,3] + LI_exact[s][1:num_days*24,4]),linewidth=3, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), #ytickfontsize=14, xtickfontsize=18,
	               legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), lc =:black, margin=5Plots.mm)
	# savefig("$dir/plots/V_Qfilter12_sum_exact.pdf")
end

control_solt = ILC_power_solt[s] + LI_exact_solt[s][:, 1] + LI_exact_solt[s][:, 2] + LI_exact_solt[s][:, 3] + LI_exact_solt[s][:, 4]
mismatch_solt = - data_solt[s] + control_solt

begin
	plot(data_solt[s], labels="Data", linecolor=:lightblue,
	fill=(0, 0, :lightblue),
	# alpha=1,
	#linestyle=:dot,
	linewidth=1,
	legend=false,
	ytickfontsize=14, xtickfontsize=18, legendfontsize=10,
	xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), margin=5Plots.mm,
	# ribbon=mismatch_solt
	)
	plot!(#control_solt,
	mismatch_solt,
	#fill=(0, 0, :lightblue),
	xticks=(0:3600*24:l_day*num_days, string.(0:num_days)),
	linecolor=:black, linestyle=:dot,
	#linewidth=1
	)
	# savefig("$dir/plots/V_Qfilter12_ILC_mismatch_secondwise_mismatch.pdf")
end

s = 1
begin
	psum_1 = plot()
	ILC_power_hourly_mean_sum =  vcat(ILC_power[s][:,:,1]'...) .+ vcat(ILC_power[s][:,:,2]'...) .+ vcat(ILC_power[s][:,:,3]'...) .+ vcat(ILC_power[s][:,:,4]'...)
	plot!(0:num_days*l_day, t -> (dd(t)[1] .- dd(t)[2]), alpha=0.7,linewidth=3, linestyle=:dot)
	plot!(1:3600:24*num_days*3600, (LI_exact[s][1:num_days*24,1] + LI_exact[s][1:num_days*24,2] + LI_exact[s][1:num_days*24,3] + LI_exact[s][1:num_days*24,4]),linewidth=3, linestyle=:dash)
	plot!(1:3600:num_days*24*3600,  ILC_power_hourly_mean_sum[1:num_days*24], xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)), #ytickfontsize=14, xtickfontsize=18,
	               legend=false, legendfontsize=10, linewidth=3,xaxis=("Days",font(10)),  yaxis=("Normed power",font(10)), lc =:black, margin=5Plots.mm)
	savefig("$dir/plots/V_noQfilter_sum_exact.png")
end

control_solt = ILC_power_solt[s] + LI_exact_solt[s][:, 1] + LI_exact_solt[s][:, 2] + LI_exact_solt[s][:, 3] + LI_exact_solt[s][:, 4]
mismatch_solt = - data_solt[s] + control_solt

begin
	plot(data_solt[s], labels="Data", linecolor=:lightblue,
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
	# savefig("$dir/plots/V_noQfilter_ILC_mismatch_secondwise_mismatch.pdf")
end
