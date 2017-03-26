/**
 * @file mgsa-core.h
 */
#ifndef __MGSA_CORE_H
#define __MGSA_CORE_H

#include "mt.h"

/**
 * Represents a prior for parameter. Very simple at the moment.
 */
struct parameter_prior
{
	/** @brief is this prior uniform continuous? (otherwise it is uniform discrete) */
	int uniform_continuous;
	double uniform_continous_lower;
	double uniform_continous_upper;

	/** @brief discrete values */
	double *values;

	/** @brief length of values */
	int number_of_states;
};

struct summary
{
	/** The parameter for which this is a summary */
	struct parameter_prior *pdsc;

	/** number of values */
	int num_of_discrete_values;

	/** The breaks, length as specified as above */
	double *breaks;

	/** The actual value counts with same length as breaks */
	int *values;

	/**
	 * Maps indices of the values field of pdsc to the corresponding
	 * index in this values. Invalid, if the described parameter
	 * is continuous.
	 */
	int *dmap;
};

int init_summary_for_breaks(struct summary *sum, int number_of_breaks);

/**
 * Represents a basic sample. Only useful in conjunction with a parameter_prior structure.
 * Can be assigned.
 */
struct prior_sample
{
	union
	{
		int discrete_index;
		double continuous_value;
	} u;
};

struct context
{
	/** @brief Number of sets */
	int number_of_sets;

	/** @brief Sizes of sets */
	int *sizes_of_sets;

	/** @brief The definitions of the sets */
	int **sets;

	/** @brief Array of size "number_of_sets" indicating the state of the sets (active or not active) */
	int *sets_active;

	/** @brief number of active sets */
	int number_of_inactive_sets;

	/**
	 *  @brief contains indices of sets that are inactive and active in form of a partition.
	 *
	 *  The pivot element is given by number_of_inactive_sets.
	 **/
	int *set_partition;

	/** @brief array that reflects the position of elements in set_inactive_list */
	int *position_of_set_in_partition;

	/** @brief Number of observable */
	int number_of_observable;

	/** @brief Array indicating whether a observable is active or not */
	int *observable;

	/** @brief Array, indicating how much sets lead to the activation of a hidden entry */
	int *hidden_count;

	/** @brief description of the prior for alpha */
	struct parameter_prior *alpha_prior;

	/** @brief prior for beta */
	struct parameter_prior *beta_prior;

	/** @brief prior for p */
	struct parameter_prior *p_prior;

	/* first digit is the observation state, second the hidden state */
	int n00;
	int n01;
	int n10;
	int n11;

	/* The actual values of the parameters */
	struct prior_sample alpha;
	struct prior_sample beta;
	struct prior_sample p;

	/* Proposal related, used to undo proposal */
	int proposal_toggle;
	int proposal_s1;
	int proposal_s2;

	int proposal_which_parameter_changed; /**! @brief defines which parameter (0-2) has been changed */
	struct prior_sample old_alpha; /* In theory, we would need only one old parameter sample (at the moment we can only change one parameter) */
	struct prior_sample old_beta;
	struct prior_sample old_p;

	/* The remainder is summary related and updated in record_activity() */
	uint64_t nsamples; /**! @brief number of samples collected */
	uint64_t *sets_activity_count;
	struct summary *alpha_summary;
	struct summary *beta_summary;
	struct summary *p_summary;

	/** @brief represents the max log likelihood encountered so far */
	double max_score;

	/** @brief in which step max_score was recorded */
	int64_t max_score_step;

	double max_score_alpha;
	double max_score_beta;
	double max_score_p;
	int *max_score_sets_active;
	int max_score_sets_active_length;
};

/** @brief MGSA MCMC parameters */
struct mcmc_params
{
	/** @brief number of steps in MCMC */
	int64_t nsteps;
	/** @brief number of burn-in MCMC steps */
	int64_t nsteps_burnin;
	/** @brief nsteps_thin sample collecting period */
	int nsteps_thin;
	/** @brief the frequency of set active state flipping Gibbs step */
	double flip_freq;
};

int init_context(struct context *cn, int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo);

void add_set(struct context *cn, int to_add);
void remove_set(struct context *cn, int to_remove);
double get_score(struct context *cn);

/* This in essence replicates the summary stuff in context. TODO: unify this */
struct result
{
	uint64_t nsamples;

	double *marg_set_activity;
	struct summary *alpha_summary;
	struct summary *beta_summary;
	struct summary *p_summary;

	double max_score;
	double max_score_alpha;
	double max_score_beta;
	double max_score_p;
	int *max_score_sets_active;
	int max_score_sets_active_length;
};

/**
 * The work horse, the single MCMC run.
 *
 * @param sets pointer to the sets. Sets a made of observable entities.
 * @param sizes_of_sets specifies the length of each set
 * @param number_of_sets number of sets (length of sets)
 * @param n the number of observable entities.
 * @param o indices of entities which are "on" (0 based).
 * @param lo length of o
 * @param params MGSA MCMC params
 * @param alpha
 * @param beta
 * @param p
 * @param mt pre-seeded structure used for random number generation.
 */
struct result do_mgsa_mcmc(int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo,
		struct mcmc_params *params,
		struct summary *alpha_summary, struct summary *beta_summary, struct summary *p_summary,
		struct mt19937p *mt);


#endif
