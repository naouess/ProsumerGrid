# Simulating demand using real load data

using JLD2, FileIO, GraphIO, CSV, DataFrames
using Distributed
using Interpolations

# here comes the broadcast
# https://docs.julialang.org/en/v1/stdlib/Distributed/index.html#Distributed.@everywhere

@everywhere begin
	dir = @__DIR__
	#include("$dir/exp_base.jl")
	include("$dir/src/experiments.jl")
#	include("$dir/input_data/demand_curves.jl")
	include("$dir/src/network_dynamics.jl")
end

@everywhere begin
		using DifferentialEquations
		using Distributions
		using LightGraphs
		using LinearAlgebra
		using Random
		using DSP
		using ToeplitzMatrices
		Random.seed!(42) # 42 random number ?
end

begin
	N = 4
	num_days =  35
	batch_size = 1
end

@everywhere begin
	freq_threshold = 0.2
	phase_filter = 1:N
	freq_filter = N+1:2N
	control_filter = 2N+1:3N
	energy_filter = 3N+1:4N
	energy_abs_filter = 4N+1:5N
end


############################################

# Defining general parameters
@everywhere begin
	l_day = 3600*24 # DemCurve.l_day
	l_hour = 3600 # DemCurve.l_hour
	l_minute = 60 # DemCurve.l_minute
	low_layer_control = experiments.LeakyIntegratorPars(M_inv=[1/5.; 1/4.8; 1/4.1; 1/4.8],kP= [400.; 110.; 100.; 200.],T_inv=[1/0.04; 1/0.045; 1/0.047; 1/0.043],kI=[0.05; 0.004; 0.05; 0.001]) # different for each node, change array
    kappa = 1. / l_hour
end

############################################
# this should only run on one process
############################################

# # Full graph for N=4 and degree 3 graph otherwise, change last 3 to 1 for N=2
_graph_lst = []
for i in 1:1
	push!(_graph_lst, random_regular_graph(iseven(3N) ? N : (N-1), 3)) # change last "3" to 1 for N=2
end
@everywhere graph_lst = $_graph_lst

####################################################
#  Preparing the demand arrays based on csv input
####################################################

struct demand_amp_var
	demand
end

function (dav::demand_amp_var)(t)
	index = Int(floor(t / (24*3600)))
	dav.demand[index + 1,:]
end

dem_data = CSV.read("$dir/profil.csv")

dem_data_week = CSV.read("$dir/profil_week.csv")

week_winter = t->dem_data_week[!,:Winterwoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
week_summer = t->dem_data_week[!,:Sommerwoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
week_between = t->dem_data_week[!,:Uebergangswoche][Int(floor(mod(t,24*3600*7) / 900)+1)]
week_G1 = t->dem_data_week[!,:WinterwocheG1][Int(floor(mod(t,24*3600*7) / 900)+1)]
week_G4 = t->dem_data_week[!,:WinterwocheG4][Int(floor(mod(t,24*3600*7) / 900)+1)]

week = t->vcat(week_winter(t), week_G1(t), week_G4(t), (week_winter(t) + week_G1(t) + week_G4(t))./3)


periodic_demand = t ->  week(t) ./100
samples = 24*60 #24*4
# Fluctuating demand! :)
inter = interpolate([10. ./100 * randn(N) for i in 1:(num_days * samples + 1)], BSpline(Linear()))
residual_demand = t -> inter(1. + t / (24*3600) * samples) # 1. + is needed to avoid trying to access out of range

#################################################
#            SIMULATION             #
# Here defining all parameters for the simulation
#################################################

# Number of ILC-Node
vc1 = 1:N # ilc_nodes (here: without communication)
# To define which nodes communicate
cover1 = Dict([v => [] for v in vc1]) # ilc_cover

# Initializing all diff. variables with zeros
u = [zeros(1000,1);1;zeros(1000,1)];

# Parameters for ILC controller
fc = 1/6;
a = digitalfilter(Lowpass(fc),Butterworth(2));
Q1 = filtfilt(a,u);#Markov Parameter

# This parameter is calculated out of the others and passed as an
# input to compound_pars later on
Q = Toeplitz(Q1[1001:1001+24-1],Q1[1001:1001+24-1]);

# Defining all parameters ! :)
_compound_pars = experiments.compound_pars(N, low_layer_control, kappa, vc1, cover1, Q)

_compound_pars.hl.daily_background_power .= 0
_compound_pars.hl.current_background_power .= 0
_compound_pars.hl.mismatch_yesterday .= 0.
_compound_pars.periodic_demand  = periodic_demand
_compound_pars.residual_demand = residual_demand
_compound_pars.graph = graph_lst[1]

# ? Why is this line needed ?
@everywhere compound_pars = $_compound_pars

# Plotting the demand
dd = t->((periodic_demand(t) .+ residual_demand(t)))
plot(0:7*l_day, t -> dd(t)[1],ytickfontsize=14,
               xtickfontsize=18, margin=5Plots.mm,
    		   legendfontsize=12, linewidth=3,xticks = (0:3600*24:num_days*24*3600, string.(0:num_days)),xaxis=("days [c]",font(14)), yaxis=("normed demand",font(14)), legend=nothing)
#title!("Demand for one week in winter (household)")
savefig("$dir/plots/real_demand_winter_week.png")

# Now putting all things together and defining the ODE
@everywhere begin
	factor = 0
	# ic means initial values :)
	ic = factor .* ones(compound_pars.D * compound_pars.N) # ? What is the parameter D ?
	tspan = (0., num_days * l_day)
	ode_tl1 = ODEProblem(network_dynamics.ACtoymodel!, ic, tspan, compound_pars,
	callback=CallbackSet(PeriodicCallback(network_dynamics.HourlyUpdate(), l_hour),
						 PeriodicCallback(network_dynamics.DailyUpdate_X, l_day)))
	# ! Callback helps retrieve current solver status for later analyses
end

@time sol1 = solve(ode_tl1, Rodas4())

# Solving 24 times per day, for number of days
hourly_energy = zeros(24*num_days+1,N)
for i=1:24*num_days+1
	# 4 Nodes
	hourly_energy[i,1] = sol1((i-1)*3600)[energy_filter[1]]
	hourly_energy[i,2] = sol1((i-1)*3600)[energy_filter[2]]
	hourly_energy[i,3] = sol1((i-1)*3600)[energy_filter[3]]
	hourly_energy[i,4] = sol1((i-1)*3600)[energy_filter[4]]
end
plot(hourly_energy)

# ILC_Power is zero for day 1 (day 0, index 1)
ILC_power = zeros(num_days+2,24,N)
# ILC_Power for day 2:
ILC_power[2,:,1] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,1])
ILC_power[2,:,2] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,2])
ILC_power[2,:,3] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,3])
ILC_power[2,:,4] = Q*(zeros(24,1) +  kappa*hourly_energy[1:24,4])
norm_energy_d = zeros(num_days,N)
norm_energy_d[1,1] = norm(hourly_energy[1:24,1])
norm_energy_d[1,2] = norm(hourly_energy[1:24,2])
norm_energy_d[1,3] = norm(hourly_energy[1:24,3])
norm_energy_d[1,4] = norm(hourly_energy[1:24,4])
mean_energy_d = zeros(num_days,N)
mean_energy_d[1,1] = mean(hourly_energy[1:24,1])
mean_energy_d[1,2] = mean(hourly_energy[1:24,2])
mean_energy_d[1,3] = mean(hourly_energy[1:24,3])
mean_energy_d[1,4] = mean(hourly_energy[1:24,4])


for i=2:num_days
	ILC_power[i+1,:,1] = Q*(ILC_power[i,:,1] +  kappa*hourly_energy[(i-1)*24+1:i*24,1])
	ILC_power[i+1,:,2] = Q*(ILC_power[i,:,2] +  kappa*hourly_energy[(i-1)*24+1:i*24,2])
	ILC_power[i+1,:,3] = Q*(ILC_power[i,:,3] +  kappa*hourly_energy[(i-1)*24+1:i*24,3])
	ILC_power[i+1,:,4] = Q*(ILC_power[i,:,4] +  kappa*hourly_energy[(i-1)*24+1:i*24,4])
	norm_energy_d[i,1] = norm(hourly_energy[(i-1)*24+1:i*24,1])
	norm_energy_d[i,2] = norm(hourly_energy[(i-1)*24+1:i*24,2])
	norm_energy_d[i,3] = norm(hourly_energy[(i-1)*24+1:i*24,3])
	norm_energy_d[i,4] = norm(hourly_energy[(i-1)*24+1:i*24,4])
	mean_energy_d[i,1] = mean(hourly_energy[(i-1)*24+1:i*24,1])
	mean_energy_d[i,2] = mean(hourly_energy[(i-1)*24+1:i*24,2])
	mean_energy_d[i,3] = mean(hourly_energy[(i-1)*24+1:i*24,3])
	mean_energy_d[i,4] = mean(hourly_energy[(i-1)*24+1:i*24,4])
end

# auswerten:
ILC_power_agg = maximum(mean(ILC_power.^2,dims=3),dims=2)
ILC_power_agg = mean(ILC_power,dims=2)
ILC_power_agg_norm = norm(ILC_power)
ILC_power_hourly = vcat(ILC_power[:,:,1]'...)
load_amp_hourly = [maximum.(dd(t)) for t in 1:3600:3600*24*num_days]

## Plots
# see original code. Will do this in a later stage
