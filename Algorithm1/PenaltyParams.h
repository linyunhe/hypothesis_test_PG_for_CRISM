#ifndef __H_PENALTY_PARAMS
#define __H_PENALTY_PARAMS

struct PenaltyParams {
	double beta_spat, delta_spat, beta_spec, delta_spec;
	bool spat_penalty_flag, spec_penalty_flag;

	int max_iter; // number of iters to apply penalty regularization
	int spectral_span;
};

#endif /* __H_PENALTY_PARAMS */