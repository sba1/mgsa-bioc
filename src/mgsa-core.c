/**
 * @file mgsa-core.c
 *
 * @author Sebastian Bauer
 *
 * This implements the core of the MGSA algorithm
 */

/*#define DEBUG*/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

#include "mgsa-core.h"
#include "mt.h"

/* Enable debugging */
//#define DEBUG

/* Define if file should be compiled as standalone (mainly for testing purposes) */
/* #define STANDALONE */

#define MIN(a,b) ((a)<(b)?(a):(b))

#ifdef STANDALONE
void *R_alloc(size_t n, int size)
{
	return malloc(n*size);
}
#else

#include <R.h>

/**
 * Thread-safe variant of R_alloc().
 *
 * @return memory allocated by R_alloc().
 * @note do not allocate memory too often.
 */
void *ts_R_alloc(size_t n, int size)
{
	void *mem;

#pragma omp critical
	{
		mem = R_alloc(n,size);
	}
	return mem;
}

/** Redirect R_alloc. I know, this is ugly... */
#define R_alloc ts_R_alloc
#endif

/*************************************************************/

/**
 * Samples from the given parameter prior.
 *
 * @param sample holds the sample
 * @param prior from which prior to sample
 * @param mt source for randomness
 *
 * @returns an opaque sample
 */
static void parameter_prior_sample(struct prior_sample *sample, struct parameter_prior *prior, struct mt19937p *mt)
{
	double rnd = genrand(mt);

	if (prior->uniform_continuous)
	{
		sample->u.continuous_value = prior->uniform_continous_lower + rnd * (prior->uniform_continous_upper - prior->uniform_continous_lower);
	} else
	{
		sample->u.discrete_index = rnd * prior->number_of_states;
		if ( sample->u.discrete_index >= prior->number_of_states )
			sample->u.discrete_index = prior->number_of_states - 1;
	}
}

/**
 * Returns the real value represented by the sample.
 *
 * @param sample
 * @param prior
 * @return
 */
static double get_prior_sample_value(struct prior_sample *sample, struct parameter_prior *prior)
{
	if (prior->uniform_continuous)
		return sample->u.continuous_value;
	return	prior->values[sample->u.discrete_index];
}


/*************************************************************/

/**
 * Allocates memory in sum to hold number_of_breaks breaks. Element
 * values is initialized while breaks is left uninitialized.
 *
 * @param sum
 * @param number_of_breaks
 * @return
 */
int init_summary_for_breaks(struct summary *sum, int number_of_breaks)
{
	sum->num_of_discrete_values = number_of_breaks;
	if (!(sum->values = R_alloc(number_of_breaks,sizeof(sum->values[0]))))
		return 0;
	memset(sum->values,0,number_of_breaks * sizeof(sum->values[0]));
	if (!(sum->breaks = R_alloc(number_of_breaks,sizeof(sum->breaks[0]))))
		return 0;
	return 1;
}

/**
 * Duplicates a given summary. The pdsc is not duplicated.
 *
 * @param sum
 * @return
 */
static struct summary *duplicate_summary(struct summary *sum)
{
	struct summary *new_sum = R_alloc(1,sizeof(*sum));
	int num_of_discrete_values = sum->num_of_discrete_values;

	if (!new_sum)
		return NULL;

	*new_sum = *sum;

	if (sum->breaks)
	{
		if (!(new_sum->breaks = (double*)R_alloc(1,sizeof(sum->breaks[0])*num_of_discrete_values)))
			return NULL;
		memcpy(new_sum->breaks,sum->breaks,sizeof(sum->breaks[0])*num_of_discrete_values);
	}

	if (sum->values)
	{
		if (!(new_sum->values = (int*)R_alloc(1,sizeof(sum->values[0])*num_of_discrete_values)))
			return NULL;
		memcpy(new_sum->values,sum->values,sizeof(sum->values[0])*num_of_discrete_values);
	}

	if (sum->dmap)
	{
		if (!(new_sum->dmap = R_alloc(num_of_discrete_values,sizeof(sum->dmap[0]))))
			return NULL;
		memcpy(new_sum->dmap,sum->dmap,sizeof(sum->dmap[0])*num_of_discrete_values);
	}
	return new_sum;
}

/**
 * Adds a new sample for the given summary.
 *
 * @param sum the summary container
 * @param value the actual value
 * @param param the description
 */
void add_to_summary(struct summary *sum, struct prior_sample *value)
{
	int slot;

	if (sum->pdsc->uniform_continuous)
	{
		double val = get_prior_sample_value(value,sum->pdsc);

		val -= sum->pdsc->uniform_continous_lower;
		/* Remember that the last slot is reserved */
		slot = val * (sum->num_of_discrete_values -1) / (sum->pdsc->uniform_continous_upper - sum->pdsc->uniform_continous_lower);
		if (slot < 0) slot=0;
		if (slot >= sum->num_of_discrete_values)
			slot = sum->num_of_discrete_values - 1;
	} else
	{
		slot = sum->dmap[value->u.discrete_index];
	}
	sum->values[slot]++;
}

/**
 * Initialize the context.
 *
 * @param context defines the context which is to be initialized
 * @param sets
 * @param sizes_of_sets
 * @param number_of_sets
 * @param n number of items
 * @param o
 * @param lo number of items that are observed as true (<=n)
 * @return
 */
int init_context(struct context *cn, int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo)
{
	int i;
#ifdef DEBUG
	printf("init_context(number_of_sets=%d,n=%d,lo=%d)\n",number_of_sets,n,lo);
#endif

	cn->number_of_sets = number_of_sets;
	cn->sets = sets;
	cn->sizes_of_sets = sizes_of_sets;
	cn->number_of_observable = n;

	/* Do the lazy stuff, i.e., allocate memory and intialize with meaningful defaults */
	if (!(cn->sets_active = (int*)R_alloc(number_of_sets,sizeof(cn->sets_active[0]))))
		goto bailout;
	memset(cn->sets_active,0,number_of_sets * sizeof(cn->sets_active[0]));
	if (!(cn->set_partition = (int*)R_alloc(number_of_sets,sizeof(cn->set_partition[0]))))
		goto bailout;
	if (!(cn->position_of_set_in_partition = (int*)R_alloc(number_of_sets,sizeof(cn->position_of_set_in_partition[0]))))
		goto bailout;
	for (i=0;i<number_of_sets;i++)
	{
		cn->set_partition[i] = i;
		cn->position_of_set_in_partition[i] = i;
	}
	cn->number_of_inactive_sets = number_of_sets;

	if (!(cn->hidden_count = (int*)R_alloc(n,sizeof(cn->hidden_count[0]))))
		goto bailout;
	memset(cn->hidden_count,0,n * sizeof(cn->hidden_count[0]));

	if (!(cn->observable = (int*)R_alloc(n,sizeof(cn->observable[0]))))
		goto bailout;
	memset(cn->observable,0,n * sizeof(cn->observable[0]));
	for (i=0;i<lo;i++)
		cn->observable[o[i]] = 1;

	if (!(cn->max_score_sets_active = (int*)R_alloc(number_of_sets,sizeof(cn->max_score_sets_active[0]))))
		goto bailout;

	/* Summary related, default values, some will be overwritten later */
	cn->nsamples = 0;
	if (!(cn->sets_activity_count = (uint64_t *)R_alloc(number_of_sets,sizeof(cn->sets_activity_count[0]))))
		goto bailout;
	memset(cn->sets_activity_count,0,number_of_sets * sizeof(cn->sets_activity_count[0]));
#if 0
	if (!(cn->alpha_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
	if (!(cn->beta_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
	if (!(cn->p_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
#endif
	/* Initially, no set is active, hence all observations that are true are false positives... */
	cn->n10 = lo;
	/* ...while the rest are true negatives */
	cn->n00 = n - lo;
	cn->n01 = cn->n11 = 0;

	cn->max_score = -DBL_MAX;
	cn->max_score_sets_active_length = 0;

	return 1;

	bailout:
	return 0;
}

/**
 * @brief A hidden member has been newly activated.
 *
 * @param cn
 * @param member
 */
static void hidden_member_activated(struct context *cn, int member)
{
	if (cn->observable[member])
	{
		/* observation is true positive rather than false positive
		 * (first digit represents the observation, the second the hidden state) */
		cn->n11++;
		cn->n10--;
	} else
	{
		/* observation is false negative rather than true negative */
		cn->n01++;
		cn->n00--;
	}
}

/**
 * @brief A hidden member has been newly deactivated.
 *
 * @param cn
 * @param member
 */
static void hidden_member_deactivated(struct context *cn, int member)
{
	if (cn->observable[member])
	{
		/* just think about it... */
		cn->n11--;
		cn->n10++;
	} else
	{
		cn->n01--;
		cn->n00++;
	}

}

/**
 * @brief add a new set.
 *
 * @param cn
 * @param to_add
 */
void add_set(struct context *cn, int to_add)
{
	int i; /* the usual dummy */

	if (cn->sets_active[to_add]) return;
	cn->sets_active[to_add] = 1;

	/* Go through all associations of the term and increase the activation count */
	for (i=0;i<cn->sizes_of_sets[to_add];i++)
	{
		int member = cn->sets[to_add][i];
		if (!cn->hidden_count[member])
		{
			/* A not yet active member is about to be activated */
			hidden_member_activated(cn,member);
		}
		cn->hidden_count[member]++;
	}

	/* Move the added set from the 0 partition to the 1 partition (it essentially becomes
	 * new first element of the 1 element, while the last 0 element gets its original position) */
	cn->number_of_inactive_sets--;
	if (cn->number_of_inactive_sets != 0)
	{
		int pos = cn->position_of_set_in_partition[to_add];
		int e0 = cn->set_partition[cn->number_of_inactive_sets];
		/* Move last element in the partition to left */
		cn->set_partition[pos] = e0;
		cn->position_of_set_in_partition[e0] = pos;
		/* Let be the newly added term the first in the partition */
		cn->set_partition[cn->number_of_inactive_sets] = to_add;
		cn->position_of_set_in_partition[to_add] = cn->number_of_inactive_sets;
	}

#ifdef DEBUG_VERBOSE
	{
		printf("Add of %d produced following sets to be inactive:\n",to_add);
		if (cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=0;i<cn->number_of_inactive_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}
		printf("Add of %d produced following sets to be active:\n",to_add);
		if (cn->number_of_sets - cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}

	}
#endif
}

/**
 * @brief remove a set.
 *
 * @param cn
 * @param to_remove
 */
void remove_set(struct context *cn, int to_remove)
{
	int i;

	if (!cn->sets_active[to_remove]) return;
	cn->sets_active[to_remove] = 0;

	/* Go through all associations of the term and decrease the activation count */
	for (i=0;i<cn->sizes_of_sets[to_remove];i++)
	{
		int member = cn->sets[to_remove][i];
		if (cn->hidden_count[member] == 1)
		{
			/* A active member is about to be deactivated  */
			hidden_member_deactivated(cn,member);
		}
		cn->hidden_count[member]--;
	}

	/* Converse of above. Here the removed set, which belonged to the 1 partition,
	 * is moved at the end of the 0 partition while the element at that place is
	 * pushed to the original position of the removed element. */
	if (cn->number_of_inactive_sets != (cn->number_of_sets - 1))
	{
		int pos = cn->position_of_set_in_partition[to_remove];
		int e1 = cn->set_partition[cn->number_of_inactive_sets];
		cn->set_partition[pos] = e1;
		cn->position_of_set_in_partition[e1] = pos;
		cn->set_partition[cn->number_of_inactive_sets] = to_remove;
		cn->position_of_set_in_partition[to_remove] = cn->number_of_inactive_sets;
	}
	cn->number_of_inactive_sets++;

#ifdef DEBUG
#ifdef DEBUG_VERBOSE
	{
		printf("Remove of %d produced following sets to be inactive:\n",to_remove);
		if (cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=0;i<cn->number_of_inactive_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}
		printf("Remove of %d produced following sets to be active:\n",to_remove);
		if (cn->number_of_sets - cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}

	}
#endif
#endif

}

/**
 * @brief Switches the state of the given set.
 *
 * @param cn defines the context on which to operate
 * @param to_switch defines the set which should be switched.
 */
static void toggle_state(struct context *cn, int to_switch)
{
	int new_state;

	new_state = !cn->sets_active[to_switch];

	if (new_state)
		add_set(cn,to_switch);
	else
		remove_set(cn,to_switch);
}

/**
 * Returns the current alpha.
 *
 * @param cn
 * @return
 */
static double get_alpha(struct context *cn)
{
	return get_prior_sample_value(&cn->alpha,cn->alpha_prior);
}

/**
 * Returns the current beta.
 *
 * @param cn
 * @return
 */
static double get_beta(struct context *cn)
{
	return get_prior_sample_value(&cn->beta,cn->beta_prior);
}

/**
 * Returns the current p.
 *
 * @param cn
 * @return
 */
static double get_p(struct context *cn)
{
	return get_prior_sample_value(&cn->p,cn->p_prior);
}

/**
 * Returns the size of the on/off neighborhood.
 * @param cn
 * @return
 */
static uint64_t get_neighborhood_size(struct context *cn)
{
	/* First part accounts for the toggles, second part for the exchanges */
	return cn->number_of_sets + cn->number_of_inactive_sets * (cn->number_of_sets - cn->number_of_inactive_sets);

}

/**
 * @brief return the score of the current context.
 *
 * @param cn
 * @return
 */
double get_score(struct context *cn)
{
	double alpha = get_alpha(cn);
	double beta = get_beta(cn);
	double p = get_p(cn);

	double score = log(alpha) * cn->n10 + log1p(-alpha) * cn->n00 + log1p(-beta)* cn->n11 + log(beta)*cn->n01;

	/* apply prior */
	score += log(p) * (cn->number_of_sets - cn->number_of_inactive_sets) + log1p(-p) * cn->number_of_inactive_sets;
	return score;
}

/**
 * @brief Proposes a new state which can be undone via undo_proposal()
 * @param cn
 * @param step current MCMC step number
 */
static void propose_state(struct context *cn, struct mcmc_params *params, struct mt19937p *mt, int64_t step)
{
	uint64_t possibilities = get_neighborhood_size(cn);

	cn->proposal_toggle = -1;
	cn->proposal_s1 = -1;
	cn->proposal_s2 = -1;

	if (genrand(mt) <= params->flip_freq)
	{
		/* toggle inactive/active states */
		uint32_t proposal = (double)(genrand(mt) * possibilities);

		if (proposal < cn->number_of_sets)
		{
			/* on/off for a single set */
			cn->proposal_toggle = proposal;
			toggle_state(cn,proposal);
		}	else
		{
			/* on/off for a pair of sets */
			int active_term_pos;
			int inactive_term_pos;

			proposal -= cn->number_of_sets;

			active_term_pos = (int)(proposal / cn->number_of_inactive_sets) +  cn->number_of_inactive_sets;
			inactive_term_pos = (int)(proposal % cn->number_of_inactive_sets);

			cn->proposal_s1 = cn->set_partition[active_term_pos];
			cn->proposal_s2 = cn->set_partition[inactive_term_pos];

			toggle_state(cn,cn->proposal_s1);
			toggle_state(cn,cn->proposal_s2);
		}
	} else
	{
		/* FIXME: If a parameter is fixed it still is considered here (this wastes steps) */
		double which_param = genrand(mt);
		if (which_param < (1.0/3.0))
		{
			cn->old_alpha = cn->alpha;
			parameter_prior_sample(&cn->alpha,cn->alpha_prior, mt);
			cn->proposal_which_parameter_changed = 0;
		} else if (which_param < (2.0/3.0))
		{
			cn->old_beta = cn->beta;
			parameter_prior_sample(&cn->beta,cn->beta_prior, mt);
			cn->proposal_which_parameter_changed = 1;
		} else
		{
			cn->old_p = cn->p;
			parameter_prior_sample(&cn->p,cn->p_prior, mt);
			cn->proposal_which_parameter_changed = 2;
		}
	}
}

/**
 * Undoes last proposal of propse_state().
 *
 * @param cn
 */
static void undo_proposal(struct context *cn)
{
	if (cn->proposal_toggle != -1) toggle_state(cn,cn->proposal_toggle);
	else if (cn->proposal_s1 != -1)
	{
		toggle_state(cn,cn->proposal_s1);
		toggle_state(cn,cn->proposal_s2);
	} else
	{
		switch (cn->proposal_which_parameter_changed)
		{
		case	0: cn->alpha = cn->old_alpha; break;
		case	1: cn->beta = cn->old_beta; break;
		default: cn->p = cn->old_p; break;
		}
	}
}

/**
 * Records the current activity of the sets. Should be called on each step.
 *
 * @param cn
 * @param score specifies the current score (log)
 */
static void record_activity(struct context *cn, int64_t step, double score)
{
	int i,j;

	cn->nsamples++;

	/* Remember that sets that are active are stored in the second partition */
	for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		cn->sets_activity_count[cn->set_partition[i]]++;

	add_to_summary(cn->alpha_summary,&cn->alpha);
	add_to_summary(cn->beta_summary,&cn->beta);
	add_to_summary(cn->p_summary,&cn->p);

	if (score > cn->max_score)
	{
		/* Remember the score and the configuration */
		cn->max_score = score;
		cn->max_score_step = step;
		cn->max_score_alpha = get_alpha(cn);
		cn->max_score_beta = get_beta(cn);
		cn->max_score_p = get_p(cn);

		for (i=cn->number_of_inactive_sets,j=0;i<cn->number_of_sets;i++,j++)
			cn->max_score_sets_active[j] = cn->set_partition[i];
		cn->max_score_sets_active_length = cn->number_of_sets - cn->number_of_inactive_sets;
	}
}

/*************************************************************/

struct result do_mgsa_mcmc(int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo,
		struct mcmc_params *params,
		struct summary *alpha_summary, struct summary *beta_summary, struct summary *p_summary,
		struct mt19937p *mt, volatile int *is_interrupted)
{
	int i;
	int64_t step;
	uint64_t neighborhood_size;
	struct context cn;
	struct result res;
	double score;

#ifdef DEBUG

	for (i=0;i<number_of_sets;i++)
	{
		int j;

		printf(" %d: ", sizes_of_sets[i]);
		for (j=0;j<sizes_of_sets[i];j++)
			printf("%d ", sets[i][j]);
		printf("\n");
	}

	printf(" n=%d\n",n);

	for (i=0;i<lo;i++)
	{
		printf(" o[%d]=%d\n",i,o[i]);
	}

	printf(" p=%d p=%d p=%d\n",alpha_prior->uniform_continuous,beta_prior->uniform_continuous,p_prior->uniform_continuous);
#endif
	memset(&res,0,sizeof(res));

	if (!init_context(&cn,sets,sizes_of_sets,number_of_sets,n,o,lo))
		goto bailout;

	if (!(alpha_summary = duplicate_summary(alpha_summary)))
		goto bailout;
	if (!(beta_summary = duplicate_summary(beta_summary)))
		goto bailout;
	if (!(p_summary = duplicate_summary(p_summary)))
		goto bailout;

	cn.alpha_summary = alpha_summary;
	cn.beta_summary = beta_summary;
	cn.p_summary = p_summary;

	cn.alpha_prior = alpha_summary->pdsc;
	cn.beta_prior = beta_summary->pdsc;
	cn.p_prior = p_summary->pdsc;

	/* Sample initial state */
	parameter_prior_sample(&cn.alpha,	cn.alpha_prior, mt);
	parameter_prior_sample(&cn.beta,	cn.beta_prior, mt);
	parameter_prior_sample(&cn.p,		cn.p_prior, mt);

	/* Disabled for now as it is not the way how the method was published
	 * This should be made optional. Note however that we expect only few
	 * terms to be on, so randomly switching on/off term seems not so
	 * appropriate */
#if 0
	for (i=0;i<number_of_sets;i++) {
		if (genrand(mt) < 0.5) toggle_state(&cn,i);
	}
#endif

	score = get_score(&cn);
	neighborhood_size = get_neighborhood_size(&cn);

#ifdef DEBUG
	printf("score=%g\n",score);
#endif

	for (step=0;step<params->nsteps;step++)
	{
		double new_score;
		double accept_probability;
		double u;
		uint64_t new_neighborhood_size;

		/* Break, if we are interrupted (we cannot use R_CheckUserInterrupt() here,
		 * as this is not thread-safe */
		if (is_interrupted && *is_interrupted) break;

		propose_state(&cn, params, mt, step);
		new_score = get_score(&cn);
		new_neighborhood_size = get_neighborhood_size(&cn);

		accept_probability = exp(new_score - score) * (double)neighborhood_size / (double)new_neighborhood_size; /* last quotient is the hasting ratio */
#ifdef DEBUG
		printf("score = %g new_score=%g alpha=%g beta=%g p=%g\n",score,new_score,get_alpha(&cn),get_beta(&cn),get_p(&cn));
#endif

		u = genrand(mt);
		if (u >= accept_probability)
		{
			undo_proposal(&cn);
		} else
		{
			score = new_score;
			neighborhood_size = new_neighborhood_size;
		}

		if ( step >= params->nsteps_burnin && (step - params->nsteps_burnin) % params->nsteps_thin == 0 ) {
			record_activity(&cn, step, score);
		}
	}

#ifdef DEBUG
	{
		int i;

		for (i=0;i<cn.number_of_sets;i++)
		{
			printf(" Set %d: %g\n",i, (double)cn.sets_activity_count[i] / (double)cn.nsamples);
		}
	}
#endif

	res.nsamples = cn.nsamples;

	if (!(res.marg_set_activity = (double*)R_alloc(number_of_sets,sizeof(res.marg_set_activity[0]))))
		goto bailout;

	for (i=0;i<cn.number_of_sets;i++)
		res.marg_set_activity[i] = (double)cn.sets_activity_count[i] / (double)cn.nsamples;

	/* Note, we are using here stuff that is allocated in init_context() */
	res.alpha_summary = cn.alpha_summary;
	res.beta_summary = cn.beta_summary;
	res.p_summary = cn.p_summary;
	res.max_score = cn.max_score;
	res.max_score_alpha = cn.max_score_alpha;
	res.max_score_beta = cn.max_score_beta;
	res.max_score_p = cn.max_score_p;
	res.max_score_sets_active = cn.max_score_sets_active;
	res.max_score_sets_active_length = cn.max_score_sets_active_length;
	bailout:
	return res;
}


static void print_context(struct context *cn)
{
	printf("n00=%d n01=%d n10=%d n11=%d num_active=%d\n",cn->n00,cn->n01,cn->n10,cn->n11,cn->number_of_sets - cn->number_of_inactive_sets);
}

#ifdef STANDALONE

int main(void)
{
	int t1[] = {0,1};
	int t2[] = {1,2};
	int o[] = {0,1};

	int *sets[] = {t1,t2};
	int sizes_of_sets[] = {sizeof(t1)/sizeof(t1[0]),sizeof(t2)/sizeof(t2[0])};
	
	struct summary alpha_summary = {0};
	struct summary beta_summary = {0};
	struct summary p_summary = {0};
	struct mt19937p mt = {0};
	struct mcmc_params params = {0};

	do_mgsa_mcmc(sets, sizes_of_sets, sizeof(sets)/sizeof(sets[0]), 3, o, sizeof(o)/sizeof(o[0]), &params, &alpha_summary, &beta_summary, &p_summary, &mt);
}

#endif
