begin
	@with_kw mutable struct η
		gen
		load
		"""
		Constructor
		"""
		function η(gen, load)
				new(gen, load)
		end
	end
	# TODO: Define outer constructors for η for the cases when only gen or load is defined

	@with_kw mutable struct LI
		kp
		ki
		T_inv
		"""
		Constructor
		"""
		function LI(kp, ki, T_inv)
				new(kp, ki, T_inv)
		end
	end

	@with_kw mutable struct ILC
		kappa
		mismatch_yesterday
		daily_background_power
		current_background_power
		ilc_nodes
		ilc_cover
		Q
		"""
		Constructor
		"""
		function ILC(kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q)
				new(kappa, mismatch_yesterday, daily_background_power, current_background_power, ilc_nodes, ilc_cover, Q)
		end
	end

	@with_kw mutable struct compound_pars
		ξ
		ILC::ILC
		LI::LI
		M_inv
		"""
		Constructor
		"""
		function compound_pars(ξ, ILC, LI, M_inv)
				new(ξ, ILC, LI, M_inv)
		end
	end
end
