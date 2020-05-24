begin
	using Random # random numbers
	using LightGraphs # create network topologies
	using LinearAlgebra
	using Parameters
	using DifferentialEquations: reinit!
	using Plots # custom plot functions
	using Interpolations
	using StatsBase
	using Statistics
	using Distributions
	using DSP
	# using DSP
	# using ToeplitzMatrices
end

l_hour = 60 * 60 # in s
l_day = l_hour * 24 # in s

@with_kw mutable struct LeakyIntegratorPars
    M_inv
    kP
    T_inv
	kI
end

@with_kw mutable struct ILCPars
	kappa
	mismatch_yesterday#::Array{Float64, 2} # Wie der Algorithmus lernt
	daily_background_power#::Array{Float64, 2} # 24xN vector with the background power for each hour.
	current_background_power
	ilc_nodes
	ilc_covers
	Q
end

@doc """
    UlmoPars(N, D, ll. hl, periodic_infeed, periodic_demand, fluctuating_infeed, residual_demand, incidence, coupling)
Define a parameter struct.
"""

# investigate this structure to define mutabe struct (with @with_kw and constructor with a function etc.)
@with_kw mutable struct UlMoPars #<: DEParameters # required parameter supertype for MCBB
	N::Int # number of nodes
	D::Int # degrees of freedom of each node
    ll::LeakyIntegratorPars # low level control parameters
	hl:: ILCPars			# high level control parameters
	periodic_infeed
	periodic_demand
	fluctuating_infeed
	residual_demand
	incidence
	coupling
	graph
	"""
	Constructor
	"""
	function UlMoPars(N::Int, #D::Int, # degrees of freedom of each node
    					ll::LeakyIntegratorPars,
						hl:: ILCPars,
						periodic_infeed,
						periodic_demand,
						fluctuating_infeed,
						residual_demand,
						#incidence,
						coupling,
						graph)
			new(N, 5,
			ll,
			hl,
			periodic_infeed,
			periodic_demand,
			fluctuating_infeed,
			residual_demand,
			incidence_matrix(graph,oriented=true),
			coupling,
			graph)
	end
end

@doc """
    default_pars(N)
Setup the system with default parameters.
"""
# ? Difference between default_pars and compound_pars?
function default_pars(N)
	low_layer_control = LeakyIntegratorPars(M_inv=60.,kP=1.3,T_inv=1.,kI=0.9)
	#g = SimpleGraph(1)
	g = random_regular_graph(iseven(3N) ? N : (N-1), 3)
	#incidence = incidence_matrix(g, oriented=true)
	# Incidence_matrix is Knoten x Kanten
	coupling= 6. .* diagm(0=>ones(size(incidence_matrix(g, oriented=true),2)))
	vc = 1:N # independent_set(g, DegreeIndependentSet()) # ilc_nodes: wo ILC sitzt
	cover = [] # Dict([v => neighbors(g, v) for v in vc]) # ilc_cover: woher Knoten Infos holen: welche Knoten kommunizieren

	# All following lines are used to calculate Q for filtering:
	u = [zeros(1000,1);1;zeros(1000,1)];
	fc = 1/6;
	a = digitalfilter(Lowpass(fc),Butterworth(2));
	Q1 = filtfilt(a,u);#Markov Parameter
	Q = Toeplitz(Q1[1001:1001+24-1],Q1[1001:1001+24-1]);

	higher_layer_control = ILCPars(kappa=0.35, mismatch_yesterday=zeros(24, N), daily_background_power=zeros(24, N), current_background_power=zeros(N), ilc_nodes=vc, ilc_covers=cover, Q=Q)
	periodic_infeed = t -> zeros(N) # Periodic infeed is zero
	peak_demand = rand(N)
	periodic_demand= t -> zeros(N)#peak_demand .* abs(sin(pi * t/24.))
	fluctuating_infeed = t -> zeros(N)
	residual_demand= t -> zeros(N)

	return UlMoPars(N, low_layer_control,
							higher_layer_control,
							periodic_infeed,
							periodic_demand,
							fluctuating_infeed,
							residual_demand,
							#incidence,
							coupling,
							g)
end

function compound_pars(N, low_layer_control, kappa, ilc_nodes, ilc_covers, Q)
	higher_layer_control = ILCPars(kappa=kappa, mismatch_yesterday=zeros(24, N), daily_background_power=zeros(24, N), current_background_power=zeros(N),ilc_nodes=ilc_nodes, ilc_covers=ilc_covers, Q=Q)

	periodic_infeed = t -> zeros(N)
	periodic_demand= t -> zeros(N)
	fluctuating_infeed = t -> zeros(N)
	residual_demand = t -> zeros(N)

	# make sure N*k is even, otherwise the graph constructor fails
	g = random_regular_graph(iseven(3N) ? N : (N-1), 3)
	coupling= 6. .* diagm(0=>ones(ne(g)))

	return UlMoPars(N, low_layer_control,
							higher_layer_control,
							periodic_infeed,
							periodic_demand,
							fluctuating_infeed,
							residual_demand,
							#incidence,
							coupling,
							g)
end

#########################################################
function ACtoymodel!(du, u, p, t)
	theta = view(u, 1:p.N)
	omega = view(u, (p.N+1):(2*p.N))
	chi = view(u, (2*p.N+1):(3*p.N))

	dtheta = view(du, 1:p.N)
	domega = view(du, (p.N+1):(2*p.N))
	dchi = view(du, (2*p.N+1):(3*p.N))

	control_power_integrator = view(du,(3*p.N+1):(4*p.N))
	control_power_integrator_abs = view(du,(4*p.N+1):(5*p.N))

	power_ILC = p.hl.current_background_power #(t)
	power_LI =  chi .- p.ll.kP .* omega
	w = p.ξ - (power_ILC + power_LI)
	@assert p.ξ .* w >= 0 # ξ and w must have the same sign, according to the defined sign convention
	@assert abs(p.ξ) .- abs(w) >= 0 # ξ always bigger than w
	u_gen = (power_ILC + power_LI) * p.η

	#periodic_power = - p.periodic_demand(t) .+ p.periodic_infeed(t)
	#fluctuating_power = - p.residual_demand(t) .+ p.fluctuating_infeed(t)

	# I don't understand this equation: how is this building the difference between the phase angles:
	flows = - (p.incidence * p.coupling * sin.(p.incidence' * theta)) # p.coupling is an important factor! Translates the line capacity

	@. dtheta = omega
	@. domega = p.ll.M_inv .* (u_gen .+ flows)
    @. dchi = p.ll.T_inv .* (- omega .- p.ll.kI .* chi) # Integrate the control power used.

	# ? these following two lines are used to aggregate. But how ?
	@. control_power_integrator = power_LI
	@. control_power_integrator_abs = abs.(power_LI)
	return nothing
end

@doc """
    HourlyUpdate()
Store the integrated control power in memory.
See also [`(hu::HourlyUpdate)`](@ref).
"""
struct HourlyUpdate
	integrated_control_power_history
	HourlyUpdate() = new([]) # needed because some definition problem occured => inner constructor! :)
end

@doc """
    HourlyUpdate(integrator)
PeriodicCallback function acting on the `integrator` that is called every simulation hour (t = 1,2,3...).
"""
function (hu::HourlyUpdate)(integrator)
	# integrator stops the solver, exracts what was te last state of variables.
	# Integrator gibt den Zustand wenn Solver angehalten wird.
	# Here the current hour is extracted using the following line:
	hour = mod(round(Int, integrator.t/3600.), 24) + 1
	# Needed for indexing (Eine Stunde geht von 0 bis 1, also letzteres nehmen)
	last_hour = mod(hour-2, 24) + 1

	# Needed for indexing
	power_idx = 3*integrator.p.N+1:4*integrator.p.N
	power_abs_idx = 4*integrator.p.N+1:5*integrator.p.N

	# the following lines: letzte Stunde auf Integral der Stunde davor setzen
	integrator.p.hl.mismatch_yesterday[last_hour,:] .= integrator.u[power_idx]
	integrator.u[power_idx] .= 0.
	integrator.u[power_abs_idx] .= 0.

	integrator.p.hl.current_background_power .= integrator.p.hl.daily_background_power[hour, :]
	nothing
end

function DailyUpdate_X(integrator)
	integrator.p.hl.daily_background_power = integrator.p.hl.Q * (integrator.p.hl.daily_background_power
	+ integrator.p.hl.kappa * integrator.p.hl.mismatch_yesterday)
	nothing
end

end
